using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using CairoMakie
using HiGHS, JSON
using JSONFBCModels

model, reaction_isozymes = build_model()

for (r,rxn) in model.reactions 
    !haskey(rxn.annotations,"REACTION") && continue 
    if any(((m,s),) -> startswith(m,"P"),rxn.stoichiometry)
        println(r)
        delete!(model.reactions,r)
    end
end

#TODO change all polymer reactions for chebis

flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
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
    ("membrane", membrane_gids, 120.0)
];
ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)




# keep membrane bound same but change biomass
ac_flux = Float64[]
vols = 0.1:0.05:3.4
bounds = Vector{Tuple{Float64,Float64}}()
for biomass in vols
    model.reactions["biomass"].upper_bound = biomass
    model.reactions["biomass"].lower_bound = biomass-0.1

    ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    push!(ac_flux,ec_sol.fluxes["EX_30089"])
    push!(bounds,(ec_sol.gene_product_capacity.membrane,ec_sol.gene_product_capacity.cytosol))
end
model.reactions["biomass"].lower_bound = 0 
model.reactions["biomass"].upper_bound = 1000
capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 400.0),
    ("membrane", membrane_gids, 60.0)
];
ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)
ac_flux_2 = Float64[]
vols_2 = 0.1:0.05:1.82
bounds_2 = Vector{Tuple{Float64,Float64}}()
for biomass in vols_2
    model.reactions["biomass"].upper_bound = biomass
    model.reactions["biomass"].lower_bound = biomass-0.1

    ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    push!(ac_flux_2,ec_sol.fluxes["EX_30089"])
    push!(bounds_2,(ec_sol.gene_product_capacity.membrane,ec_sol.gene_product_capacity.cytosol))
end

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

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
)
lines!(
    ax,
    vols_2,
    abs.(round.(ac_flux_2))./vols_2,
    label = "Membrane: 60mg/gDW",
    linewidth=2.5

)
lines!(
    ax,
    vols,
    abs.(round.(ac_flux))./vols,
    label = "Membrane: 120mg/gDW",
    color = :red,
    linestyle = :dash,
    linewidth=2.5

)
xlims!(ax,(0,3.7))
axislegend(
    ax,
    position=:lt,
    labelsize = 5pt,
)
display(f)

vols[337]
ac_flux[337]

save("data/plots/acetate_exchange.png",f,px_per_unit = 1200/inch)
