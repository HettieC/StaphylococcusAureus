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
using Latexify
using ColorSchemes
flux_zero_tol = 1e-6
gene_zero_tol = 1e-6

model, reaction_isozymes = build_model()

gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)

model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_15903"].upper_bound = 0 #block glucose exchange
model.reactions["EX_16651"].upper_bound = 1000 #allow l-lactate uptake

model.reactions["EX_16480"] = CM.Reaction(;
    name="nitric oxide exchange",
    stoichiometry=Dict(
        "CHEBI:16480" => 1, #nitric oxide
    ),
    lower_bound = 0.1,
)

capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 120.0),
    ("membrane", membrane_gids, 120.0)
];

ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)

# exchanges
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

open("data/fluxes.json","w") do io
    JSON.print(io, ec_sol.fluxes)
end

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
pruned_sol.gene_product_capacity

# differentiation variables
parameter_values = Dict(
    Symbol(x) => y[iso_id].kcat_forward for (x,y) in pruned_reaction_isozymes for (iso_id,iso) in y if ec_sol.isozyme_forward_amounts[x][iso_id] > gene_zero_tol || ec_sol.isozyme_reverse_amounts[x][iso_id] > gene_zero_tol
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
for r in parameters[order][203:end]
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
    reverse(sens[flux_idx,:][order][203:end]),
    color=Makie.wong_colors()[2]
)
ylims!(ax,(0,maximum(reverse(sens[flux_idx,:][order][203:end]))*1.05))
f

save("data/plots/NO_biomass_sens.png", f, px_per_unit = 1200/inch)

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

save("data/plots/NO_resp_sens.png", f, px_per_unit = 1200/inch)

# look at whole solution, excluding permeases
biomass_id = :biomass 
biomass_idx = findfirst(x -> last(x) == biomass_id && first(x) == :fluxes, vids)
flux_idxs = findall(x -> !startswith(string(last((x))),"PERM") && first(x) == :fluxes, vids)
flux_ids = last.(vids[flux_idxs])
flux_order = sortperm(lowercase.(string.(flux_ids)))
push!(flux_order,length(flux_ids)+1) # add biomass to the end of the order
push!(flux_idxs, biomass_idx) # add biomass to the end of the flux indices
push!(flux_ids, biomass_id) # add biomass to the end of the flux ids

param_idxs = findall(x -> !startswith(string(x),"PERM"), parameters)
param_ids = parameters[param_idxs]

xtickvals = []
for i in axes(sens[flux_idxs[flux_order], param_idxs]',1)
    if abs(sum(sens[flux_idxs[flux_order], param_idxs]'[i,:])/length(flux_ids))>0.02
        push!(xtickvals,i)
    end
end
xtickvals
xticklabels = String[] 
for r in string.(parameters[param_order][xtickvals])
    if isnothing(tryparse(Int,string(r)))
         push!(xticklabels,A.reaction_name(model,string(r)))
    elseif !haskey(model.reactions[string(r)].annotations,"BiGG") || isempty(model.reactions[string(r)].annotations["BiGG"])
        push!(xticklabels,string(r)*(length(A.reaction_name(escher_model,string(r))) > 15 ? "" : ": "*A.reaction_name(escher_model,string(r))))
    else 
        push!(xticklabels,string(r)*": "*model.reactions[string(r)].annotations["BiGG"][1])
    end
end
xticklabels

sub_sens = zeros(size(sens[flux_idxs[flux_order], param_idxs]' ))
for (i,j) in zip(1:size(sens[flux_idxs[flux_order], param_idxs]',1),1:size(sens[flux_idxs[flux_order], param_idxs]',2))
    sub_sens[i,j] = sens[flux_idxs[flux_order], param_idxs]'[i,j] < 0 ? -log(abs(sens[flux_idxs[flux_order], param_idxs]'[i,j])) : log(sens[flux_idxs[flux_order], param_idxs]'[i,j])
end

f = Figure()#, backgroundcolor=:transparent)
ax = Axis(
    f[1,1],
    xlabel = "Turnover number",
    xticklabelrotation = -pi / 2,
    ylabel = "Flux",
    xticks = (xtickvals, xticklabels),
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
)
hm = heatmap!(
    ax,
    sens[flux_idxs[flux_order], param_idxs]';
    colormap = reverse(ColorSchemes.RdBu),
    colorrange = (-0.02, 0.02)
)
Colorbar(f[1, 2], hm, ticklabelsize = 5pt)
f

