using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
import COBREXA as X
using CairoMakie
using HiGHS, JSON, CSV
using JSONFBCModels, DataFrames
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
using CairoMakie
using Latexify
using ColorSchemes
flux_zero_tol = 1e-6
gene_zero_tol = 1e-6

model, reaction_isozymes = build_model()

gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)

# double the kcat of permease reactions 
for (k,v) in reaction_isozymes 
    !startswith(k,"PERM") && continue 
    for (id,iso) in v 
        reaction_isozymes[k][id].kcat_forward *= 2
        reaction_isozymes[k][id].kcat_reverse *= 2
    end
end

model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange

capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 200.0),
    ("membrane", membrane_gids, 130.0)
];

ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)

# atp-producing rxns
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:30616") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)



pruned_model, pruned_reaction_isozymes = D.prune_model(
    model,
    ec_sol.fluxes,
    ec_sol.gene_product_amounts,
    reaction_isozymes,
    ec_sol.isozyme_forward_amounts,
    ec_sol.isozyme_reverse_amounts,
    flux_zero_tol,
    gene_zero_tol,
);

pruned_sol = enzyme_constrained_flux_balance_analysis(
    pruned_model;
    reaction_isozymes = pruned_reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)

Dict(x=>[y,abs(ec_sol.fluxes[x])] for (x,y) in pruned_sol.fluxes if abs(abs(y)-abs(ec_sol.fluxes[x]))>1e-5)


# differentiation variables
parameter_values = Dict(
    Symbol(x) => y[iso_id].kcat_forward for (x,y) in pruned_reaction_isozymes for (iso_id,iso) in y if ec_sol.isozyme_forward_amounts[x][iso_id] > 1e-8 || ec_sol.isozyme_reverse_amounts[x][iso_id] > 1e-8
)

parameters_2 = Symbol.(collect(keys(parameter_values)))
rid_pid = Dict(rid => Ex(Symbol(rid)) for (rid, v) in pruned_reaction_isozymes)
rid_gcounts = Dict(rid => [v.gene_product_stoichiometry for (k, v) in d][1] for (rid, d) in pruned_reaction_isozymes)

parameter_isozymes = Dict{String,Dict{String,X.IsozymeT{Ex}}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
parameter_isozymes = Dict(
    rid => Dict(
        iso_id => X.IsozymeT{Ex}(
            gene_product_stoichiometry = iso.gene_product_stoichiometry,
            kcat_forward = rid_pid[rid],
            kcat_reverse = nothing
        )
        for (iso_id,iso) in y
    ) 
    for (rid,y) in pruned_reaction_isozymes
)

pkm = X.enzyme_constrained_flux_balance_constraints( # pruned kinetic model
    pruned_model;
    reaction_isozymes = parameter_isozymes,
    gene_product_molar_masses,
    capacity
)

pruned_solution = D.optimized_values(
    pkm,
    parameter_values;
    objective = pkm.objective.value,
    optimizer = HiGHS.Optimizer,
)

pkm_kkt, vids = D.differentiate_prepare_kkt(pkm, pkm.objective.value, parameters_2)
sens_2 = D.differentiate_solution(
    pkm_kkt,
    pruned_solution.primal_values,
    pruned_solution.equality_dual_values,
    pruned_solution.inequality_dual_values,
    parameter_values,
    scale = true, # unitless sensitivities
)


# look at biomass 
flux_id = [:biomass]
flux_idx_2 = findall(x -> last(x) == :biomass && first(x) == :fluxes, vids)

order_2 = sortperm(sens_2[flux_idx, :],dims=2)

sens_2[flux_idx_2, :][1:end]
inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=15)

# combined plot

sens_dict = Dict(string.(parameters) .=> sens[flux_idx,:]')
sens_dict_2 = Dict(string.(parameters_2) .=> sens_2[flux_idx_2,:]')
params = [string(p) for p in reverse(parameters[order])[1:20]] 
unique!(append!(params,[string(p) for p in reverse(parameters_2[order_2])[1:20]]))

data = (
    xs = params,
    ys1 = float.([haskey(sens_dict,p) ? sens_dict[p] : 0 for p in params]),
    ys2 = float.([haskey(sens_dict_2,p) ? sens_dict_2[p] : 0 for p in params]),
)
xticklabels = String[] 
for r in data.xs
    if isnothing(tryparse(Int,string(r)))
         push!(xticklabels,A.reaction_name(model,string(r)))
    elseif !haskey(model.reactions[string(r)].annotations,"BiGG") || isempty(model.reactions[string(r)].annotations["BiGG"])
        push!(xticklabels,string(r)*(length(A.reaction_name(escher_model,string(r))) > 15 ? "" : ": "*A.reaction_name(escher_model,string(r))))
    else 
        push!(xticklabels,string(r)*": "*model.reactions[string(r)].annotations["BiGG"][1])
    end
end

xticklabels
f = Figure(; size=(12cm, 9cm))#, backgroundcolor=:transparent)

ax = Axis(
    f[1,1],
    #backgroundcolor=:transparent,
    xlabel = "Enzyme",
    xticklabelrotation = -pi / 3,
    ylabel = "Sensitivity of biomass",
    xlabelsize=7pt,
    ylabelsize=7pt,
    xticklabelsize=6pt,
    yticklabelsize=6pt,
    xticks = (1:length(data.xs),xticklabels),
    ygridvisible=false,
    xgridvisible=false,
)

barplot!(
    ax,
    repeat(1:length(data.xs),2),
    vcat(data.ys1,data.ys2),
    dodge = vcat(repeat([1],length(data.xs)),repeat([2],length(data.xs))),
    color = vcat(repeat([Makie.wong_colors()[2]],length(data.xs)),repeat([Makie.wong_colors()[6]],length(data.xs))),
)
ylims!(ax,(0,0.23))
xlims!(ax, (0, length(data.xs)+1))
labels = ["60/h","120/h"]
elements = [
    PolyElement(color=Makie.wong_colors()[2]),
    PolyElement(color=Makie.wong_colors()[6]),
]
Legend(f[1,1],
elements,
labels,
labelsize=8pt,
tellheight = false,
tellwidth = false,
margin = (10, 10, 10, 10),
halign = :right, valign = :top,
)
f

save("data/plots/biomass_sens_combined.png", f, px_per_unit = 1200/inch)


[r for r in params if haskey(model.reactions[r].annotations,"Pathway") && "Glycolysis / Gluconeogenesis; " in model.reactions[r].annotations["Pathway"]]


rs = [string(p) for p in reverse(parameters_2[order_2])[1:10]]
for r in rs 
    println(r)
    println(A.reaction_name(model, r))
    for (k,v) in model.reactions[r].annotations
        println("  $k: $v")
    end
    println("")
end

for r in rs 
    if "Pentose phosphate pathway; " ∈ model.reactions[r].annotations["Pathway"]
        println(r)
        println(A.reaction_name(model, r))
        for (k,v) in model.reactions[r].annotations
            println("  $k: $v")
        end
        println("")
    end
end
