using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using CairoMakie
using HiGHS, JSON
using JSONFBCModels

model, reaction_isozymes = build_model();
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)
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
model.reactions["biomass"].objective_coefficient = 0
ac_flux = Float64[]
growth = Float64[]
vols = 0.1:0.05:copy(ec_sol.objective)+0.033
bounds = Vector{Tuple{Float64,Float64}}()
for biomass in vols
    model.reactions["biomass"].upper_bound = biomass
    model.reactions["biomass"].lower_bound = biomass-0.02

    ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity)
    
    ec_sol = optimized_values(
        ct;
        optimizer = HiGHS.Optimizer,
        objective = ct.fluxes.EX_30089.value,
        sense = Maximal
    )
    push!(ac_flux,abs(ec_sol.fluxes["EX_30089"]))
    push!(bounds,(ec_sol.gene_product_capacity.membrane,ec_sol.gene_product_capacity.cytosol))
    push!(growth,ec_sol.fluxes["biomass"])
end

model.reactions["biomass"].objective_coefficient = 1
model.reactions["biomass"].lower_bound = 0 
model.reactions["biomass"].upper_bound = 1000
capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 200.0),
    ("membrane", membrane_gids, 110.0)
];
ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)
model.reactions["biomass"].objective_coefficient = 0

ac_flux_2 = Float64[]
growth_2 = Float64[]
vols_2 = 0.1:0.05:copy(ec_sol.objective)
bounds_2 = Vector{Tuple{Float64,Float64}}()
for biomass in vols_2
    model.reactions["biomass"].upper_bound = biomass
    model.reactions["biomass"].lower_bound = biomass-0.001

    ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity)
    
    ec_sol = optimized_values(
        ct;
        optimizer = HiGHS.Optimizer,
        objective = ct.fluxes.EX_30089.value,
        sense = Maximal
    )
    push!(ac_flux_2,abs(ec_sol.fluxes["EX_30089"]))
    push!(bounds_2,(ec_sol.gene_product_capacity.membrane,ec_sol.gene_product_capacity.cytosol))
    push!(growth_2,ec_sol.fluxes["biomass"])
end

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)
colors = Makie.wong_colors()

f = Figure(; size=(10cm, 6cm))

ax = Axis(
    f[1,1];
    backgroundcolor=:transparent,
    xlabel = "Growth rate (gDW/h)",
    ylabel = "Acetate exchange rate (mMol/gDW/h)",
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    #xticks = [0,0.5,1,1.5,2,2.5],
    ygridvisible=false,
    xgridvisible=false,
);
band!(
    ax,
    growth[findfirst(((x,y),) -> x>119.999,bounds)]:0.0001:growth[end],
    0,
    maximum(ac_flux),
    color = (colors[1],0.3)
)
band!(
    ax,
    growth_2[findfirst(((x,y),) -> x>109.999,bounds_2)]:0.0001:growth_2[end],
    0,
    maximum(ac_flux),
    color = (colors[6],0.3)
)
lines!(
    ax,
    growth_2,
    ac_flux_2,
    label = "Membrane: 110mg/gDW",
    linewidth=2.5,
    color = colors[6]
)
lines!(
    ax,
    growth,
    ac_flux,
    label = "Membrane: 120mg/gDW",
    color = colors[1],
    linestyle = :dash,
    linewidth=2.5
)
xlims!(ax,(0,growth[end]))
ylims!(ax, (0,maximum(ac_flux)))
axislegend(
    ax,
    position=:lt,
    labelsize = 5pt,
)
display(f)

growth[findfirst(x -> x> 1e-3,ac_flux)]
growth_2[findfirst(x -> x> 1e-3,ac_flux_2)]

save("data/plots/acetate_min_acetate.png",f,px_per_unit = 1200/inch)
