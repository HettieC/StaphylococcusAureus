using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using CairoMakie
using HiGHS, JSON
using JSONFBCModels

model, reaction_isozymes = build_model()

model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_15903"].upper_bound = 10 #limit glucose

sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
sol.fluxes["EX_15379"]
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

open("data/fluxes.json","w") do io 
    JSON.print(io,sol.fluxes)
end

# scan growth rate across limited oxygen with fba model 
growth = Float64[]
o2_iter = 0:1:35
for o2_uptake in o2_iter
    model.reactions["EX_15379"].upper_bound = o2_uptake
    model.reactions["EX_15379"].lower_bound = o2_uptake-0.1
    sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    push!(growth,sol.objective)
end
growth

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

#atp producing
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:30616") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)

# oxygen scan with EC model 
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_15903"].upper_bound = 1000 #unlimit glucose
model.reactions["EX_15379"].lower_bound = 0 #unlimit o2
model.reactions["EX_15379"].upper_bound = 1000

capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 200.0),
    ("membrane", membrane_gids, 120.0)
];
ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
) 
ec_sol.fluxes["EX_15379"]
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)
ec_growth = Float64[]
bounds = Vector{Tuple{Float64,Float64}}()
for o2_uptake in o2_iter
    println(o2_uptake)
    model.reactions["EX_15379"].upper_bound = o2_uptake
    model.reactions["EX_15379"].lower_bound = o2_uptake-0.1
    ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )    
    isnothing(ec_sol) && break
    push!(ec_growth,ec_sol.objective)
    push!(bounds,(ec_sol.gene_product_capacity.cytosol,ec_sol.gene_product_capacity.membrane))
end
ec_growth
bounds
inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

f = Figure(; size=(10cm, 6cm))#, backgroundcolor=:transparent)

colors = Makie.wong_colors()

ax = Axis(
    f[1,1];
    backgroundcolor=:transparent,
    ylabel = "Growth rate (gDW/h)",
    xlabel = "Oxygen uptake rate (mmol/gDW/h)",
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    ygridvisible=false,
    xgridvisible=false,
    xticks = [0,10,20,o2_iter[findfirst(x -> x==maximum(growth),growth)],o2_iter[findfirst(((x,y),)->y>119.999,bounds)],30]
)
band!(
    ax,
    o2_iter[findfirst(((x,y),)->y>119.999,bounds)]:maximum(o2_iter)+0.1,
    0,
    maximum(vcat(growth,ec_growth))*1.1,
    color = (colors[6],0.3)
)
lines!(
    ax,
    o2_iter[1:length(ec_growth)],
    ec_growth,
    label = "ecFBA",
    linewidth = 2.5,
    color = colors[6]
)
lines!(
    ax,
    o2_iter,
    growth,
    label = "FBA",
    color=colors[1],
    linestyle=:dash,
    linewidth=2.5
)
vlines!(
    ax,
    o2_iter[findfirst(x -> x==maximum(growth),growth)],
    color = colors[4],
    linestyle = :dash,
    linewidth = 0.5
)
axislegend(
    ax,
    position=:cb,
    labelsize = 5pt,
)
xlims!(ax,(0,maximum(o2_iter)+0.1))
ylims!(ax,(0,maximum(vcat(growth,ec_growth))*1.1))
f
save("data/plots/o2_scan.png", f, px_per_unit = 1200/inch)

