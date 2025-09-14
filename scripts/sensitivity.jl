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

escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")
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

parameters = Symbol.(collect(keys(parameter_values)))
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

pkm_kkt, vids = D.differentiate_prepare_kkt(pkm, pkm.objective.value, parameters)
sens = D.differentiate_solution(
    pkm_kkt,
    pruned_solution.primal_values,
    pruned_solution.equality_dual_values,
    pruned_solution.inequality_dual_values,
    parameter_values,
    scale = true, # unitless sensitivities
)


# look at biomass 
flux_id = [:biomass]
flux_idx = findall(x -> last(x) == :biomass && first(x) == :fluxes, vids)

order = sortperm(sens[flux_idx, :],dims=2)

sens[flux_idx, :][1:end]
inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=12)

f = Figure(; size=(10cm, 8cm))#, backgroundcolor=:transparent)
escher_model = change_reaction_names(pruned_model)
xticklabels = String[] 
for r in parameters[order][212:end]
    if isnothing(tryparse(Int,string(r)))
         push!(xticklabels,string(r))
    elseif !haskey(model.reactions[string(r)].annotations,"BiGG") || isempty(model.reactions[string(r)].annotations["BiGG"])
        push!(xticklabels,string(r))
    else 
        push!(xticklabels,string(r)*": "*model.reactions[string(r)].annotations["BiGG"][1])
    end
end

ax = Axis(
    f[1,1],
    #backgroundcolor=:transparent,
    xlabel = "Enzyme",
    xticklabelrotation = -pi / 3,
    ylabel = "Sensitivity of biomass",
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    xticks = (1:length(xticklabels),reverse(xticklabels)),
    ygridvisible=false,
    xgridvisible=false,
)
barplot!(
    ax,
    1:length(xticklabels),
    reverse(sens[flux_idx,:][order][212:end]),
    color=Makie.wong_colors()[2]
)
ylims!(ax,(0,0.23))
f

save("data/plots/biomass_sens.png", f, px_per_unit = 1200/inch)



# look at ox phos 
subset_ids = [:cyt_bo3, :Sdh, :Mqo, :ATPS, :Ndh2, :Lqo, :Ldh, :cyt_aa3, :cyt_bd]
biomass_id = :biomass 
biomass_idx = findfirst(x -> last(x) == biomass_id && first(x) == :fluxes, vids)
flux_idxs = findall(x -> last(x) in subset_ids && first(x) == :fluxes, vids)
flux_ids = last.(vids[flux_idxs])
flux_order = sortperm(lowercase.(string.(flux_ids)))
push!(flux_order,length(flux_ids)+1) # add biomass to the end of the order
push!(flux_idxs, biomass_idx) # add biomass to the end of the flux indices
push!(flux_ids, biomass_id) # add biomass to the end of the flux ids

param_idxs = findall(x -> x in subset_ids, parameters)
param_ids = parameters[param_idxs]
param_order = sortperm(lowercase.(string.(param_ids)))

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

f = Figure(; size=(10cm, 8cm))#, backgroundcolor=:transparent)
ax = Axis(
    f[1,1],
    xlabel = L"\text{Enzyme, }p",
    xticks = (1:length(param_ids), string.(param_ids[param_order])),
    xticklabelrotation = -pi / 2,
    ylabel = L"Reaction flux sensitivity, $\frac{\partial v}{\partial p}$",
    yticks = (1:length(flux_ids), string.(flux_ids[flux_order])),
    xlabelsize=7pt,
    ylabelsize=7pt,
    xticklabelsize=6pt,
    yticklabelsize=6pt,
)
hm = heatmap!(
    ax,
    sens[flux_idxs[flux_order], param_idxs[param_order]]';
    colormap = reverse(ColorSchemes.RdBu[1:6]),
    colorrange = (0,maximum(sens[flux_idxs[flux_order], param_idxs[param_order]]))
)
Colorbar(f[1, 2], hm, ticklabelsize = 5pt)
f

save("data/plots/respiration_sens.png", f, px_per_unit = 1200/inch)



# whole solution

using Clustering

flux_idxs = findall(x -> first(x) == :fluxes, vids)
flux_ids = last.(vids[flux_idxs])
# cluster flux sensitivity into 10 clusters using K-means
R = kmeans(sens[flux_idxs,:]',20; maxiter=400,display=:iter)

a = assignments(R) # get the assignments of points to clusters

# make heatmap of sens[flux_idxs[sortperm(a)],:]' with the flux_ids ordered by their parameter_location

xtickvals = []
for i in axes(sens[flux_idxs[sortperm(a)],:]',1)
    if (sum(abs.(sens[flux_idxs[sortperm(a)],:]'[i,:]))/length(parameters))>0.06
        push!(xtickvals,i)
    end
end
xtickvals
xticklabels = [isnothing(A.reaction_name(escher_model,p)) ? A.reaction_name(model,p) : A.reaction_name(escher_model,p) for p in string.(parameters[xtickvals])]


xtickvals = findfirst(x->string(x)=="PTS_CHEBI:15903",parameters)
xticklabels = string.(parameters[xtickvals])

ytickvals = []
for i in axes(sens[flux_idxs[sortperm(a)],:]',2)
    if (sum(abs.(sens[flux_idxs[sortperm(a)],:]'[:,i]))/length(flux_ids))>0.1
        push!(ytickvals,i)
    end
end
ytickvals
yticklabels = [length(A.reaction_name(escher_model,p)) > 10 ? p : A.reaction_name(escher_model,p) for p in string.(flux_ids[ytickvals])]


f = Figure(; size=(12cm,10cm))#, backgroundcolor=:transparent)
ax = Axis(
    f[1,1],
    xlabel = L"\text{Enzyme, }p",
    xticklabelrotation = -pi / 3,
    ylabel = L"Reaction flux sensitivity, $\frac{\partial v}{\partial p}$",    
    xticks = ([xtickvals], [xticklabels]),
    xlabelsize=7pt,
    ylabelsize=7pt,
    xticklabelsize=6pt,
    yticklabelsize=6pt,
    yticksvisible = false,
    yticklabelsvisible = false,
    #yticks = (ytickvals, yticklabels),
)
hm = heatmap!(
    ax,
    sens[flux_idxs[sortperm(a)][201:end],:]';
    colormap = reverse(ColorSchemes.RdBu),
    colorrange = (-0.2,0.2)
)
Colorbar(f[1, 2], hm, ticklabelsize = 5pt)
f

save("data/plots/whole_sens.png", f, px_per_unit = 1200/inch)


