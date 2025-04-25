using DataFrames, CSV
using CairoMakie
using Statistics
using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
using COBREXA, JSONFBCModels
using HiGHS, JSON
using ConstraintTrees
import ConstraintTrees as C
using Serialization

model = build_model()

reaction_isozymes = deserialize("reaction_isozymes.serialized")

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

gene_product_molar_masses = get_gene_product_molar_mass([g for g in A.genes(model) if g!="g1"])
gene_product_molar_masses["g1"] = mean(collect(values(gene_product_molar_masses)))


total_capacity = 550.0 # kDa

# simulate overflow metabolism 
all_gids = A.genes(model)
membrane_gids = ["g1"]

membrane_capacity = 0.009 * total_capacity


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

refsim = run_ecfba_monoculture(ct, 100)

# make plots

# acetate overflow 
wt = scatter(
    refsim.mu,
    abs.(refsim.ac_flux),;
    axis = (xlabel = "Growth rate [1/h]", 
            ylabel = "Acetate flux [gDW/mmol/h]", 
            xlabelsize= 25, 
            ylabelsize = 25, 
            xticklabelsize = 20,
            yticklabelsize = 20),
)

open("data/refsim_fluxes_growth_0.775.json","w") do io 
    JSON.print(io,refsim.flux_sol[3])
end

save("wt_ac_tiny_mem.png",wt)