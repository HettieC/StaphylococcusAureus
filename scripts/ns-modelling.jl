using StaphylococcusAureus
import AbstractFBCModels as A
import ConstraintTrees as C
using COBREXA
using DataFrames, CSV
using CairoMakie
using Statistics
using HiGHS

model = build_model()

turnup_df = DataFrame(CSV.File("data/model/isozymes/turnup_output.csv"))

rxn_seq_subs_prods = DataFrame(CSV.File("data/model/isozymes/reaction_sequence_subs_prods.csv"))

df = DataFrame(CSV.File("data/model/isozymes/locus_tag_seq_ST398.csv"))
g_id_sequence = Dict(Pair.(df.first, df.second))

turnup_df = DataFrames.rename!(turnup_df,
    "kcat [s^(-1)]" => "kcat",
)
kcat_df = insertcols(rxn_seq_subs_prods, 5, :kcat => turnup_df.kcat)
# kcat mean at 22.82
kcat_dict = Dict{String,Dict{String,Float64}}()
for row in eachrow(kcat_df)
    if !ismissing(row.kcat)
        if !haskey(kcat_dict, row.reaction_id)
            # put into 1/h instead of 1/s by multiplying by 3600
            kcat_dict[row.reaction_id] = Dict(row.locus_tag => row.kcat * 3.6)
        else
            kcat_dict[row.reaction_id][row.locus_tag] = row.kcat * 3.6
        end
    end
end
avg_kcat = mean(vcat([collect(values(y)) for (x,y) in kcat_dict]...))

isozymes_stoich_df = DataFrame(CSV.File("data/model/isozymes/reaction_isozymes.csv")) 
select!(isozymes_stoich_df,:RHEA_ID, :Protein, :Stoichiometry, :Isozyme)
isozyme_stoich = Dict(Pair.(isozymes_stoich_df.Protein, isozymes_stoich_df.Stoichiometry))

reaction_isozymes = Dict{String,Dict{String,Isozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for (id, rxn) in model.reactions
    grrs = rxn.gene_association_dnf
    isnothing(grrs) && continue # skip if no grr available
    haskey(kcat_dict, "$(id)_f") || continue # skip if no kcat data available
    haskey(kcat_dict, "$(id)_r") || continue
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, id, Dict{String,Isozyme}())
        #println(id,"  ",[isozyme_stoich[g] for g in grr])
        d["isozyme_"*string(i)] = Isozyme( # each isozyme gets a unique name
            gene_product_stoichiometry = Dict(g => isozyme_stoich[g] for g in grr), 
            kcat_forward = maximum([kcat_dict["$(id)_f"][g] for g in grr]),
            kcat_reverse = maximum([kcat_dict["$(id)_r"][g] for g in grr]),
        )
    end
end


# add oxphos fake isozymes 
oxphos_reactions = ["Ndh2", "Sdh", "Mqo", "Lqo", "cyt_aa3", "cyt_bd", "cyt_bo3"]
for rid in oxphos_reactions
    grrs = A.reaction_gene_association_dnf(model, rid)
    reaction_isozymes[rid] = Dict("isoyzme_1" => Isozyme(
        gene_product_stoichiometry = Dict("g1" => 1), # assume subunit stoichiometry of 1 for all isozymes
        kcat_forward = avg_kcat,
        kcat_reverse = avg_kcat,
        )
    )
end


gene_product_molar_masses = get_gene_product_molar_mass([g for g in A.genes(model) if g!="g1"])
gene_product_molar_masses["g1"] = mean(collect(values(gene_product_molar_masses)))


total_capacity = 550.0 # kDa


ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses = gene_product_molar_masses,
    capacity = total_capacity,
    optimizer = HiGHS.Optimizer,
)
# ec_sol.gene_product_amounts
# # The total amount of required gene product mass
# ec_sol.gene_product_capacity

# open("data/ec_fluxes.json","w") do io 
#     JSON.print(io,ec_sol.fluxes)
# end

# simulate overflow metabolism 
all_gids = A.genes(model)
membrane_gids = ["g1"]

membrane_capacity = 0.01 * total_capacity


function prepare_ec_ecoli(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    all_gids,
    total_capacity,
    membrane_gids,
    membrane_capacity,
)
    ct = enzyme_constrained_flux_balance_constraints(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity = [
            ("total", all_gids, total_capacity),
            ("membrane", membrane_gids, membrane_capacity),
        ],
    )

    # apply more quirks (make a few bounds more realistic)
    ct.fluxes[:EX_15903].bound = C.Between(-1000, 1000)
    ct.fluxes_forward[:EX_15903].bound = C.Between(0, 1000)
    ct.fluxes_reverse[:EX_15903].bound = C.Between(0, 1000)
    # ct.fluxes[:EX_5dglcn_e].bound = C.EqualTo(0.0)
    # ct.fluxes[:EX_for_e].bound = C.EqualTo(0.0)
    # ct.fluxes[:EX_pyr_e].bound = C.EqualTo(0.0)
    # ct.fluxes[:EX_lac__D_e].bound = C.EqualTo(0.0)

    # inject the overexpressed proteins into the system
    # ct +=
    #     :gene_product_amounts^C.variables(; keys = [:lacY, :lacZ], bounds = C.EqualTo(0.0))
    # ct.gene_product_capacity.total.value +=
    #     ct.gene_product_amounts.lacY.value * lacY_mm +
    #     ct.gene_product_amounts.lacZ.value * lacZ_mm
    # ct.gene_product_capacity.membrane.value += ct.gene_product_amounts.lacY.value * lacY_mm

    # add the "total amount of gene stuff" objective that we want to minimize
    ct *=
        :l1_min_proteins_objective^C.Constraint(
            sum(
                C.value(v) * gene_product_molar_masses[string(k)] for
                (k, v) in ct.gene_product_amounts
            ),
            nothing,
        )

    return ct
end

function run_ecfba_monoculture(ct, length; lacZ_frac = 0.0, lacY_frac = 0.0)
    # make a copy for local modifications
    ct = C.ConstraintTree(
        ct...,
        :objective => deepcopy(ct.objective),
        :gene_product_amounts => deepcopy(ct.gene_product_amounts),
    )

    # ct.gene_product_amounts[:lacY].bound = C.EqualTo(total_capacity * lacY_frac / lacY_mm)
    # ct.gene_product_amounts[:lacZ].bound = C.EqualTo(total_capacity * lacZ_frac / lacZ_mm)

    res = screen(range(0.1, 1.2, length)) do mu
        #@info "ecfba run" mu lacY_frac lacZ_frac

        ct.objective.bound = C.EqualTo(mu)

        sol = optimized_values(
            ct;
            objective = ct.l1_min_proteins_objective.value,
            optimizer = HiGHS.Optimizer,
            settings = [silence],
            sense = Minimal,
        )

        isnothing(sol) && return nothing

        return (;
            mu,
            total_mass = sol.gene_product_capacity.total,
            membrane_mass = sol.gene_product_capacity.membrane,
            ac_flux = sol.fluxes.EX_30089,
            glc_flux = sol.fluxes.EX_15903,
            o2_flux = sol.fluxes.EX_15379,
            ATP_flux = sol.fluxes.ATPS,
            flux_sol = sol.fluxes,
        )
    end
    return Tables.columntable(filter(!isnothing, res))
end


ct = prepare_ec_ecoli(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    all_gids,
    total_capacity,
    membrane_gids,
    membrane_capacity,
)

refsim = run_ecfba_monoculture(ct, 10)

wt = scatter(
    refsim.mu,
    abs.(refsim.ac_flux);
    axis = (xlabel = "Growth rate [1/h]", 
            ylabel = "Acetate flux [mmol/gDW/h]", 
            xlabelsize=20, 
            ylabelsize = 20, 
            xticklabelsize = 20,
            yticklabelsize = 20),
)

open("data/refsim_fluxes_growth_0.775.json","w") do io 
    JSON.print(io,refsim.flux_sol[3])
end

save("wt_ac.png",wt)
