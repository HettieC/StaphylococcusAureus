using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using DataFrames, CSV
using CairoMakie
using Statistics
using HiGHS, XLSX, DataFramesMeta
using ElementaryFluxModes, JSON
import DifferentiableMetabolism as D
using JSONFBCModels

# add this to transporters.csv: Permease,glucose,CHEBI:15903,SAPIG2309,1
model, reaction_isozymes = build_model()
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")
model.reactions["EX_16236"].lower_bound = 0 #block ethanol exchange
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange

# add fake isozymes 
oxphos_reactions = ["Ndh2", "Sdh", "Mqo", "Lqo", "cyt_aa3", "cyt_bd", "cyt_bo3","Ldh"]
for rid in A.reactions(model)
    haskey(reaction_isozymes,rid) && continue
    if rid ∈ oxphos_reactions || !isnothing(tryparse(Int,rid))
        println(rid)
        grrs = A.reaction_gene_association_dnf(model, rid)
        reaction_isozymes[rid] = Dict("isoyzme_$i" => Isozyme(
            gene_product_stoichiometry=Dict(g => 1 for g in grr), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward=65,
            kcat_reverse=65,
        )
        for (i,grr) in enumerate(model.reactions[rid].gene_association_dnf)
        )
    end
end


gene_product_molar_masses = get_gene_product_molar_mass([g for g in A.genes(model) if g != "g1"])
gene_product_molar_masses["g1"] = mean(collect(values(gene_product_molar_masses)))

cytosol_gids = String[] 
for (r,isos) in reaction_isozymes 
    if !isnothing(tryparse(Int,r))
        append!(cytosol_gids,[g for (x,y) in isos for (g,s) in y.gene_product_stoichiometry])
    end
end
membrane_gids = [g for g in A.genes(model) if g ∉ cytosol_gids]

df1 = DataFrame(CSV.File("data/databases/uniprot/s_aureus.csv"))
df2 = DataFrame(CSV.File("data/databases/uniprot/idmapping_2025_05_19.csv"))
rename!(df2,:Entry => :uniprot_accesion)

df = innerjoin(df1,df2,on=:uniprot_accesion)

subcellular_location = Dict{String,Vector{String}}()
for ln in eachrow(df)        
    if ismissing(ln.Subcellular_location_CC)
        subcellular_location[ln.gene] = ["uncategorised"]
    else
        subs = string.(split(split(ln.Subcellular_location_CC,"SUBCELLULAR LOCATION: ")[2],"; "))
        subcellular_location[ln.gene] = rstrip.(first.(split.(subs," {")),'.')
    end
end
append!(membrane_gids, [x for (x,y) in subcellular_location if x ∈ A.genes(model) && any(z -> occursin("membrane",lowercase(z)),y) ])

capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 300.0),
    ("membrane", membrane_gids, 10.0)
]

# model.reactions["EX_30089"].objective_coefficient = 0
# model.reactions["biomass"].objective_coefficient = 1
# model.reactions["biomass"].stoichiometry["CHEBI:30089"] = -0.01

ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)

open("data/fluxes.json","w") do io 
    JSON.print(io,ec_sol.fluxes)
end
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)


ac_flux = Float64[]
biomass = Float64[]
vols = 1:2:100
for vol in vols
    capacity = [
        ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 300.0),
        ("membrane", membrane_gids, vol)
    ]
    ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    push!(ac_flux,ec_sol.fluxes["EX_30089"])
    push!(biomass,ec_sol.fluxes["biomass"])
end


using CairoMakie
fig = Figure(;size = (700,500))
ax = Axis(
    fig[1,1],
    xlabel = "Membrane bound, mg/gDW",
    ylabel = "Flux"
)
lines!(
    ax,
    vols,
    abs.(ac_flux),
    label = "Acetate exchange"
)
lines!(
    ax,
    vols,
    abs.(biomass),
    label = "Biomass reaction"
)
fig[1,2] = Legend(fig,ax)
fig



# keep membrane bound same but change biomass
ac_flux = Float64[]
membrane_conc = Float64[]
vols = 0:0.05:5
for biomass in vols
    model.reactions["biomass"].upper_bound = biomass
    ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    push!(ac_flux,ec_sol.fluxes["EX_30089"])
    push!(membrane_conc,ec_sol.gene_product_capacity.membrane)
end
fig = Figure(;size = (700,500))
ax = Axis(
    fig[1,1],
    xlabel = "Growth rate",
    ylabel = "Acetate exchange"
)
lines!(
    ax,
    vols,
    abs.(ac_flux),
    label = "Acetate exchange"
)
lines!(
    ax,
    vols,
    abs.(membrane_conc),
    label = "Membrane conc"
)
fig[1,2] = Legend(fig,ax)
fig


model.reactions["biomass"].upper_bound = 2.9
ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
)

open("data/fluxes.json","w") do io 
    JSON.print(io,ec_sol.fluxes)
end


C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 100    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

# prune the model 
flux_zero_tol = 1e-6 # these bounds make a real difference!
gene_zero_tol = 1e-6
pruned_model, pruned_reaction_isozymes = D.prune_model(
    model,
    ec_sol.fluxes,
    ec_sol.gene_product_amounts,
    reaction_isozymes,
    ec_sol.isozyme_forward_amounts,
    ec_sol.isozyme_reverse_amounts,
    flux_zero_tol,
    gene_zero_tol,
)
pruned_solution = enzyme_constrained_flux_balance_analysis(
    pruned_model;
    reaction_isozymes = pruned_reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)

N = A.stoichiometry(pruned_model)
atpm_idx = findfirst(x -> x == "ATPM", A.reactions(pruned_model))
biomass_idx = findfirst(x -> x == "biomass", A.reactions(pruned_model))
fixed_fluxes = [atpm_idx, biomass_idx]
flux_values = [
    pruned_solution.fluxes["ATPM"],
    pruned_solution.fluxes["biomass"],
]


OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)
OFM_dicts = [
    Dict(A.reactions(pruned_model) .=> OFMs[:, 1]),
    Dict(A.reactions(pruned_model) .=> OFMs[:, 2]),
]
Dict([x,model.reactions[x].name] => [y,OFM_dicts[2][x]] for (x,y) in OFM_dicts[1] if abs(y - OFM_dicts[2][x])>1)


