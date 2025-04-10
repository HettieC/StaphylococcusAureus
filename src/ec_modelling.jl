using COBREXA
using DataFrames, CSV
using CairoMakie

turnup_df = DataFrame(CSV.File("data/model/isozymes/turnup_output.csv"))

rxn_seq_subs_prods = DataFrame(CSV.File("data/model/isozymes/reaction_sequence_subs_prods.csv"))

df = DataFrame(CSV.File("data/model/isozymes/locus_tag_seq_ST398.csv"))
g_id_sequence = Dict(Pair.(df.first, df.second))

turnup_df = DataFrames.rename!(turnup_df,
    "kcat [s^(-1)]" => "kcat",
)
kcat_df = insertcols(rxn_seq_subs_prods, 5, :kcat => turnup_df.kcat)
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
            gene_product_stoichiometry = Dict(g => isozyme_stoich[g] for g in grr), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = maximum([kcat_dict["$(id)_f"][g] for g in grr]),
            kcat_reverse = maximum([kcat_dict["$(id)_r"][g] for g in grr]),
        )
    end
end

function get_gene_product_molar_mass(gids)
    AA_mass = Dict(
        'A' => 89,
        'R' => 174,
        'N' => 132,
        'D' => 133,
        'B' => 133,
        'C' => 121,
        'Q' => 146,
        'E' => 147,
        'Z' => 147,
        'G' => 75,
        'H' => 155,
        'I' => 131,
        'L' => 131,
        'K' => 146,
        'M' => 149,
        'F' => 165,
        'P' => 115,
        'S' => 105,
        'T' => 119,
        'W' => 204,
        'Y' => 181,
        'V' => 117,
    )
    seqs = Dict{String,String}()
    open("data/sequence.txt","r") do io
        prev = ""
        for ln in eachline(io)
            if startswith(ln,'>')
                seqs[ln[2:end]] = ""
                prev = ln[2:end]
            else 
                seqs[prev] = seqs[prev]*ln 
            end
        end
    end

    gene_product_molar_mass = Dict{String,Float64}()
    for gid in gids 
        gene_product_molar_mass[gid] = sum([AA_mass[aa] for aa in seqs[gid]]) / 1000 # convert to milligrams
    end

    println(110*length(seqs["SAPIG1968"]))
    return gene_product_molar_mass     
    
end

gene_product_molar_masses = get_gene_product_molar_mass(A.genes(model))

total_capacity = 550.0 # kDa

ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses = gene_product_molar_masses,
    capacity = total_capacity,
    optimizer = HiGHS.Optimizer,
)
ec_sol.gene_product_amounts
# The total amount of required gene product mass
ec_sol.gene_product_capacity

open("data/ec_fluxes.json","w") do io 
    JSON.print(io,ec_sol.fluxes)
end

# simulate overflow metabolism 
all_gids = A.genes(model)




membrane_capacity = 0.20 * total_capacity

function prepare_ec_ecoli(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    all_gids,
    total_capacity,
    # membrane_gids,
    # membrane_capacity,
)
    ct = enzyme_constrained_flux_balance_constraints(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity = [
            ("total", all_gids, total_capacity),
            #("membrane", membrane_gids, membrane_capacity),
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

    res = screen(range(0.1, 0.24, length)) do mu
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
            # membrane_mass = sol.gene_product_capacity.membrane,
            ac_flux = sol.fluxes.EX_30089,
            glc_flux = sol.fluxes.EX_15903,
            o2_flux = sol.fluxes.EX_15379,
            ATP_flux = sol.fluxes.ATPS,
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
    # membrane_gids,
    # membrane_capacity,
)

refsim = run_ecfba_monoculture(ct, 20)

wt = scatter(
    refsim.mu,
    abs.(refsim.o2_flux);
    axis = (xlabel = "Growth rate [1/h]", ylabel = "Acetate flux [mmol/gDW/h]", xlabelsize=20, ylabelsize = 20),
)

open("data/refsim_fluxes.json","w") do io 
    JSON.print(io,ec_sol.fluxes)
end