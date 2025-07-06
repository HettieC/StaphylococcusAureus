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
using Latexify
using LaTeXStrings

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
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 400.0),
    ("membrane", membrane_gids, 110.0)
];


ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)
# ec_sol.gene_product_capacity

# open("data/fluxes.json","w") do io 
#     JSON.print(io,ec_sol.fluxes)
# end
# C.pretty(
#     C.ifilter_leaves(ec_sol.fluxes) do ix, x
#         abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
#     end; 
#     format_label = x -> A.reaction_name(model, string(last(x))),
# )
# # atp-producing rxns
# C.pretty(
#     C.ifilter_leaves(ec_sol.fluxes) do ix, x
#         abs(x) > 1e-5 && 
#             haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:30616") && 
#             ((model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] > 0 && x > 1e-5) || 
#             (model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] < 0 && x < -1e-5))
#     end; 
#     format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
# )


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

# diff_flux = Dict(string(x)=>[y-abs(ec_sol.fluxes[x])] for (x,y) in pruned_sol.fluxes if abs(abs(y)-abs(ec_sol.fluxes[x]))>1e-5)

#pruned_sol.gene_product_capacity

# sum([pruned_sol.gene_product_amounts[g]*gene_product_molar_masses[g] for g in membrane_gids if haskey(pruned_sol.gene_product_amounts,g)])
# sum([e*gene_product_molar_masses[string(g)] for (g,e) in pruned_sol.gene_product_amounts if string(g) ∉ membrane_gids])

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

OFM_dicts[1]["ATPS"]
OFM_dicts[2]["ATPS"]

Dict(x => [y,OFM_dicts[2][x]] for (x,y) in OFM_dicts[1] if startswith(x,"EX"))

ofm_df = DataFrame(Exchange_metabolite=String[],OFM_1=Float64[],OFM_2=Float64[])
for (x,y) in OFM_dicts[1]
    if startswith(x,"EX")
        push!(
            ofm_df,
            [
                A.metabolite_name(model,"CHEBI:"*string(split(x,"_")[2])),
                round(y,sigdigits=4),
                round(OFM_dicts[2][x],sigdigits=4)
            ])
    end
end
ofm_df
sort!(ofm_df,:OFM_1)
latexify(ofm_df; env = :table, booktabs = true, latex = false) |> print

# df of big difference reactions 
ofm_df = DataFrame(Reaction=String[],Name=String[],OFM_1=Float64[],OFM_2=Float64[])
for (x,y) in OFM_dicts[1]
    if abs(y - OFM_dicts[2][x])/y > 1 || abs(y - OFM_dicts[2][x])/OFM_dicts[2][x] > 0.75
        push!(
            ofm_df,
            [
                x,
                isnothing(A.reaction_name(model,x)) ? A.reaction_name(escher_model,x) : A.reaction_name(model,x),
                round(y,sigdigits=4),
                round(OFM_dicts[2][x],sigdigits=4)
            ]
        )
    end
end
ofm_df
sort!(ofm_df,:OFM_1)
latexify(ofm_df; env = :table, booktabs = true, latex = false) |> print



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

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)
f = Figure(; size = (18cm,9cm))

data = (
    x=1:length(parameters),
    height1=[c[1] > 0 ? c[1] : c[2] for c in eachcol(scaled_sens[:, order])],
    height2=abs.([c[1] < 0 ? c[1] : c[2] for c in eachcol(scaled_sens[:, order])]),
    grp1=[c[1] > 0 ? 1 : 2 for c in eachcol(scaled_sens[:, order])],
    grp2=[c[1] < 0 ? 1 : 2 for c in eachcol(scaled_sens[:, order])]
    #grp1 = [x == "cytosol" ? 1 : 2 for x in grp],
    #grp2 = [x == "membrane" ? 2 : 1 for x in grp],
)
ax1 = Axis(
    f[1, 1], 
    xreversed = true, 
    xautolimitmargin = (0, 0.1),
    xlabelsize=9pt,
    ylabelsize=6pt,
    xticklabelsize=8pt,
    xgridvisible=false,
    ygridvisible=false,
    yticksvisible = false,
    xscale = log10,
)
ax2 = Axis(
    f[1, 2], 
    alignmode = Mixed(left = Makie.Protrusion(0)), 
    xautolimitmargin = (0, 0.1),
    xlabelsize=9pt,
    ylabelsize=6pt,
    xticklabelsize=8pt,
    yticklabelsize=8pt,
    xgridvisible=false,
    ygridvisible=false,
    xscale = log10,
    yticksvisible = false,
    yticks = ([findlast(x -> x == 2, data.grp1)/2,length(filter(g->g==1,data.grp1))/2+findlast(x -> x == 2, data.grp1)],[L"\textbf{Membrane   }\;", L"\textbf{Cytosol   }\;"]),
    #yticklabelrotationion = π/2
)
hideydecorations!(ax1, grid = false)
linkyaxes!(ax1, ax2)
colgap!(f.layout, 0)
barplot!(ax1, data.x, data.height1, direction = :x, color = colors[data.grp1], label = "OFM 1: Respiratory")
barplot!(ax2, data.x, data.height2, direction = :x, color = colors[data.grp2], label = "OFM 2: Fermentative")
xlims!(ax1, [10,8e-6])
xlims!(ax2, [2e-6,20])
labels = [L"\textbf{Respiratory OFM}", L"\textbf{Fermentative OFM}"]
elements = [MarkerElement(marker=:hline,color=colors[2],markersize=12), MarkerElement(marker=:hline,color=colors[1],markersize=12)]
Legend(
    f[1, 1],
    elements,
    labels,
    tellheight=false,
    tellwidth=false,
    labelsize = 8pt,
    halign = :left,
    valign = :center,
    margin = (16,0,0,0),
    framevisible = false
)
bracket!(ax1, 8e-6, 0, 8e-6, findlast(x -> x == 2, data.grp1), style=:curly, orientation=:up,linewidth=1, width = 10)
bracket!(ax1, 8e-6, findlast(x -> x == 2, data.grp1), 8e-6, length(data.grp1), style=:curly, orientation=:up,linewidth=1, width = 10)
Label(f[2,:],L"\textbf{OFM Sensitivity: }\frac{p}{\lambda}\frac{\partial \lambda}{\partial p}")
f

save("data/plots/ofm.png", f, px_per_unit = 1200/inch)
