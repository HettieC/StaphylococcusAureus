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

flux_zero_tol = 1e-10
gene_zero_tol = 1e-10

model, reaction_isozymes = build_model();

gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)

escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange

model.reactions["EX_15903"].upper_bound = 1000

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
open("data/fluxes.json","w") do io 
    JSON.print(io,ec_sol.fluxes)
end
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
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

Dict(x => [y,ec_sol.fluxes[x]] for (x,y) in pruned_sol.fluxes if abs(ec_sol.fluxes[x]) - y >1e-5)


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
                round(y,sigdigits=3),
                round(OFM_dicts[2][x],sigdigits=3)
            ])
    end
end
ofm_df
sort!(ofm_df,:OFM_1)
latexify(ofm_df; env = :table, booktabs = true, latex = false) |> print

# df of big difference reactions 
ofm_df = DataFrame(Reaction=String[],Name=String[],OFM_1=Float64[],OFM_2=Float64[])
for (x,y) in OFM_dicts[1]
    if y > 1.2*OFM_dicts[2][x] || OFM_dicts[2][x] > 1.2*y
        name = isnothing(A.reaction_name(model,x)) ? A.reaction_name(escher_model,x) : A.reaction_name(model,x) 
        isnothing(name) && println(x)
        isnothing(name) && continue
        push!(
            ofm_df,
            [
                x,
                isnothing(A.reaction_name(model,x)) ? A.reaction_name(escher_model,x) : A.reaction_name(model,x),
                round(y,sigdigits=3),
                round(OFM_dicts[2][x],sigdigits=3)
            ]
        )
    end
end

ofm_df
sort!(ofm_df,:OFM_1)
latexify(ofm_df; env = :table, booktabs = true, latex = false) |> print

# OFM fluxes for escher:
open("data/ofm1_fluxes.json","w") do io 
    JSON.print(io,OFM_dicts[1])
end
open("data/ofm2_fluxes.json","w") do io 
    JSON.print(io,OFM_dicts[2])
end
escher_model = change_reaction_names(pruned_model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/pruned_escher_model.json")



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
colors = Makie.wong_colors()[3:4]

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
)
# left hand side
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
    xticks = ([1,0.01,0.0001],[L"-10^0",L"-10^{-2}",L"-10^{-4}"]),
)
#right hand side
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
    xticks = ([0.0001,0.01,1],[L"10^{-4}",L"10^{-2}",L"10^0",]),
    xscale = log10,
    yticksvisible = false,
    yticks = ([findlast(x -> x == 2, data.grp1)/2,length(filter(g->g==1,data.grp1))/2+findlast(x -> x == 2, data.grp1)],[L"\textbf{Membrane   }\;",L"\textbf{Cytosol   }\;"]),
    #yticklabelrotationion = π/2
)
display(f)
hideydecorations!(ax1, grid = false)
linkyaxes!(ax1, ax2)
colgap!(f.layout, 0)
barplot!(ax1, data.x, data.height2, direction = :x, color = colors[data.grp2], label = "OFM 1: Fermentative")
barplot!(ax2, data.x, data.height1, direction = :x, color = colors[data.grp1], label = "OFM 2: Respiratory")
xlims!(ax1, [8,2e-5])
xlims!(ax2, [2e-5,20])
f 
labels = [L"\textbf{Fermentative OFM}",L"\textbf{Respiratory OFM}",]
elements = [MarkerElement(marker=:hline,color=colors[1],markersize=12), MarkerElement(marker=:hline,color=colors[2],markersize=12)]
Legend(
    f[1, 1],
    elements,
    labels,
    tellheight=false,
    tellwidth=false,
    labelsize = 8pt,
    position = (0,0.3),
    halign = :left,
    valign = :bottom,
    margin = (16,0,50,0),
    framevisible = false
)
bracket!(ax1, 2e-5, 0, 2e-5, findlast(x -> x == 2, data.grp1), style=:curly, orientation=:up,linewidth=1, width = 10)
bracket!(ax1, 2e-5, findlast(x -> x == 2, data.grp1), 2e-5, length(data.grp1), style=:curly, orientation=:up,linewidth=1, width = 10)
Label(f[2,:],L"\textbf{OFM Sensitivity: }\frac{p}{\lambda}\frac{\partial \lambda}{\partial p}")
display(f)

save("data/plots/ofm.png", f, px_per_unit = 1200/inch)

findfirst(x->x=="ATPS",string.(parameters[order]))

df = DataFrame(Parameter=String[],Parameter_name=String[],OFM_1=Float64[],OFM_2=Float64[])
for (i,p) in enumerate(parameters)
    push!(df, [string(p), isnothing(A.reaction_name(model,string(p))) ? "" : A.reaction_name(model,string(p))  ,scaled_sens[1, order[i]], scaled_sens[2, order[i]]])
end
df
sort!(df,:OFM_1)

df2 = filter(row -> row.Parameter ∈ ofm_df.Reaction,df)



#############
## check flux_ids[sortperm(a)] are used in the same quantitiy by the two OFMs
to_check = string.(flux_ids[sortperm(a)][201:end])
Dict(x => [y,OFM_dicts[2][x]] for (x,y) in OFM_dicts[1] if x in to_check && abs(y-OFM_dicts[2][x]) < 1e-10)

most_controlling = parameters[order][vcat(1:5,222:226...)]
[A.reaction_name(escher_model,string(r)) for r in most_controlling]

dic = Dict(r => [i,A.reaction_name(escher_model,string(r)),haskey(model.reactions[string(r)].annotations,"Pathway") ? model.reactions[string(r)].annotations["Pathway"] : "" ] for (i,r) in enumerate(string.(most_controlling)))

findfirst(x->x==1,data.grp1)
data.height2[29:35]

rs = string.(parameters[order][29:35])
dic = Dict(r => [A.reaction_name(escher_model,r),haskey(model.reactions[r].annotations,"Pathway") ? model.reactions[r].annotations["Pathway"] : "" ] for r in rs)

sens = Dict(string(x) => i ∈ 29:35 ? 1 : 0 for (i,x) in enumerate(string.(parameters[order])) )
for r in A.reactions(model)
    if !haskey(sens,r)
        sens[r] = 0
    end
end
open("data/sens.json","w") do io 
    JSON.print(io,Dict(string.(parameters[order])[i] => data.height1[i] for i in 1:length(data.height1)))
end

open("data/sens_bold.json","w") do io 
    JSON.print(io,sens)
end
