using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
import COBREXA as X
using ElementaryFluxModes
using CairoMakie
using HiGHS, JSON, CSV
using JSONFBCModels, DataFrames
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
using CairoMakie

flux_zero_tol = 1e-6
gene_zero_tol = 1e-6

model, reaction_isozymes = build_model()

gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)

escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")
model.reactions["EX_16236"].lower_bound = 0 #block ethanol exchange
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_16651"].lower_bound = 0 #block (s)-lactate exchange
model.reactions["EX_16004"].lower_bound = 0 #block (r)-lactate exchange
model.reactions["EX_15740"].lower_bound = 0 #block formate exchange
model.reactions["EX_15378"].lower_bound = 0 #block H+ exchange


capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 300.0),
    ("membrane", membrane_gids, 120.0)
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




#### calculate efms 
# calculate EFMs
N = A.stoichiometry(pruned_model)

# atpm, biomass
atpm_idx = findfirst(x -> x == "ATPM", A.reactions(pruned_model))
biomass_idx = findfirst(x -> x == "biomass", A.reactions(pruned_model))
fixed_fluxes = [atpm_idx, biomass_idx]
flux_values = [pruned_sol.fluxes["ATPM"], pruned_sol.fluxes["biomass"]]

OFMs = get_ofms(Matrix(N), fixed_fluxes, flux_values)

OFM_dicts = [
    Dict(A.reactions(pruned_model) .=> OFMs[:,1]),
    Dict(A.reactions(pruned_model) .=> OFMs[:,2])
]
# scale to have one unit flux through biomass
OFM_dicts = [
    Dict(x=>y/OFM_dicts[1]["biomass"] for (x,y) in OFM_dicts[1]),
    Dict(x=>y/OFM_dicts[2]["biomass"] for (x,y) in OFM_dicts[2])
]

OFM_dicts[1]["EX_30089"]/OFM_dicts[1]["EX_15903"]
OFM_dicts[2]["EX_30089"]/OFM_dicts[2]["EX_15903"]

OFM_dicts[1]["EX_15379"] #o2
OFM_dicts[2]["EX_15379"]

Dict(x => [y,OFM_dicts[2][x]] for (x,y) in OFM_dicts[1] if startswith(x,"EX"))


## calculate lambda
nonzero_1 = [x for (x, y) in OFM_dicts[1] if OFM_dicts[2][x] == 0]
nonzero_2 = [x for (x, y) in OFM_dicts[2] if OFM_dicts[1][x] == 0]
M = [
    OFM_dicts[1][nonzero_1[1]] OFM_dicts[2][nonzero_1[1]];
    OFM_dicts[1][nonzero_2[1]] OFM_dicts[2][nonzero_2[1]]
]

v = [
    pruned_sol.fluxes[nonzero_1[1]];
    pruned_sol.fluxes[nonzero_2[1]]
]

# rounding causes issues
lambda = inv(M) * v

parameter_values = Dict(
    Symbol(x) => y[iso_id].kcat_forward for (x,y) in pruned_reaction_isozymes for (iso_id,iso) in y if ec_sol.isozyme_forward_amounts[x][iso_id] > 1e-8 || ec_sol.isozyme_reverse_amounts[x][iso_id] > 1e-8
)
parameters = Ex.(collect(keys(parameter_values)))
rid_pid = Dict(rid => Ex(Symbol(rid)) for (rid, v) in pruned_reaction_isozymes)
rid_gcounts = Dict(rid => [v.gene_product_stoichiometry for (k, v) in d][1] for (rid, d) in pruned_reaction_isozymes)

sens_efm = differentiate_efm(OFM_dicts, parameters, rid_pid, parameter_values, rid_gcounts, capacity, gene_product_molar_masses, HiGHS.Optimizer)

param_vals = collect(values(parameter_values))
scaled_sens = Matrix(undef, size(sens_efm, 1), size(sens_efm, 2))
for (i, col) in enumerate(eachcol(sens_efm))
    scaled_sens[:, i] = param_vals[i] .* col ./ lambda
end

order = sortperm(scaled_sens[1, :])
colors = Makie.wong_colors()[5:6]

# ensure plot is correct size
inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=1)

fig = Figure(; size=(10cm, 7cm))#, backgroundcolor=:transparent)
data = (
    x=1:length(parameters),
    height=[c[1] > 0 ? c[1] : c[2] for c in eachcol(scaled_sens[:, order])],
    grp=[c[1] > 0 ? 1 : 2 for c in eachcol(scaled_sens[:, order])],
)
ax = Axis(
    fig[1, 1],
    ylabel="OFM control coefficient",
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    xgridvisible=false,
    ygridvisible=false,
    #xlabel = "Parameter",
    xticks=([length([c for c in eachcol(scaled_sens) if c[1]<0])/2, size(scaled_sens,2)-length([c for c in eachcol(scaled_sens) if c[1]>0])/2], ["Cytosolic\nEnzymes", "Membrane\nEnzymes"]),
    yscale=log10
)
barplot!(
    ax,
    data.x,
    data.height,
    color=colors[data.grp]
)
labels = ["More fermentative OFM", "More respiratory OFM"]
elements = [MarkerElement(marker=:rect,color=colors[2],markersize=10), MarkerElement(marker=:rect,color=colors[1],markersize=10)]
Legend(
    fig[1, 1],
    elements,
    labels,
    tellheight=false,
    tellwidth=false,
    labelsize = 5pt,
    halign = :center,
    valign = :top,
    margin = (20,0,0,0),
    framevisible = false
)
bracket!(ax, 0, 4e-7, findlast(x -> x == 2, data.grp), 4e-7, style=:curly, orientation=:down,linewidth=0.5, width = 10)
bracket!(ax, findfirst(x -> x == 1, data.grp), 4e-7, length(data.grp), 4e-7, style=:curly, orientation=:down, linewidth=0.5,width=10)
fig

save("data/plots/ofm.png", fig, px_per_unit = 1200/inch)

### location of parameter 
param_loc = Dict(rid => any(g -> g ∈ membrane_gids, collect(keys(first(iso).second.gene_product_stoichiometry))) ? "membrane" : "cytosol" for (rid,iso) in pruned_reaction_isozymes)

[param_loc[string(x)] for x in parameters[order]]

filter(((x,y),) -> y=="membrane",param_loc)

[r for (r,v) in pruned_sol.fluxes if isnothing(tryparse(Int,string(r)))]

