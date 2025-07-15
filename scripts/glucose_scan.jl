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
model.reactions["EX_15379"].upper_bound = 10 #limit oxygen
model.reactions["EX_15379"].lower_bound = 0
model.reactions["EX_15378"].upper_bound = 0
model.reactions["EX_15378"].lower_bound = 0 #block H+ exchange


sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
sol.fluxes["EX_15903"]

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

# scan growth rate across limited oxygen with fba model 
growth = Float64[]
glc_rate = Float64[]
glc_iter = 0:10:560
for glc in glc_iter
    model.reactions["EX_15903"].upper_bound = glc
    model.reactions["EX_15903"].lower_bound = glc-0.1
    sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    push!(growth,sol.objective)
    push!(glc_rate,glc/sol.objective)
end
growth
glc_rate

# oxygen scan with EC model 
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)
model.reactions["EX_16236"].lower_bound = 0 #block ethanol exchange
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_16651"].lower_bound = 0 #block (s)-lactate exchange
model.reactions["EX_16004"].lower_bound = 0 #block (r)-lactate exchange
model.reactions["EX_15903"].upper_bound = 1000 #unlimit glucose
model.reactions["EX_15903"].lower_bound = 0 
model.reactions["EX_15379"].upper_bound = 1000 #unlimit o2
model.reactions["EX_15378"].lower_bound = 0 #block H+ exchange


capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 400.0),
    ("membrane", membrane_gids, 120.0)
];

ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
) 
ec_sol.fluxes["EX_15903"]

ec_growth = Float64[]
bounds = Vector{Tuple{Float64,Float64}}()
for glc in 1:2:40
    println(glc)
    model.reactions["EX_15903"].upper_bound = glc
    model.reactions["EX_15903"].lower_bound = glc-0.1
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
    glc/ec_sol.objective > maximum(glc_rate) && break
end
ec_growth
bounds
inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

f = Figure(; size=(10cm, 6cm))#, backgroundcolor=:transparent)

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
    #xticks = [0,5,10,15]
)
lines!(
    ax,
    glc_iter[1:length(ec_growth)]./ec_growth,
    ec_growth,
    label = "ecFBA",
    linewidth = 2.5
)
lines!(
    ax,
    glc_rate,
    growth,
    label = "FBA",
    color=:red,
    linestyle=:dash,
    linewidth=2.5
)
axislegend(
    ax,
    position=:cb,
    labelsize = 5pt,
)
xlims!(ax,(0,maximum(o2_rate)+0.1))
ylims!(ax,(0,maximum(vcat(growth,ec_growth))*1.1))
f
save("data/plots/glc_scan_gDW.png", f, px_per_unit = 1200/inch)
