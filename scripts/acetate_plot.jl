using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using CairoMakie
using HiGHS, JSON
using JSONFBCModels

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
    ("membrane", membrane_gids, 120.0)
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

C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:30616") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)


using CairoMakie
# keep membrane bound same but change biomass
ac_flux = Float64[]
membrane_conc = Float64[]
growth = []
vols = 0.1:0.1:ec_sol.objective
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
    push!(growth,ec_sol.objective) 
    push!(membrane_conc,ec_sol.gene_product_capacity.membrane)
end
ac_flux
growth
membrane_conc


inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

f = Figure(; size=(10cm, 6cm))#, backgroundcolor=:transparent)

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
    vols,
    abs.(ac_flux)./growth,
    label = "Acetate exchange"
)
axislegend(
    ax,
    position=:lt,
    labelsize = 5pt,
)
f


save("plots/acetate_exchange.png",f,px_per_unit = 1200/inch)


model.reactions["biomass"].upper_bound = 1000

# iterate over membrane bounds
ac_flux = Float64[]
membrane_capacity = Float64[]
vols = 20:10:300
for membrane_bound in vols
    capacity = [
        ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 200.0),
        ("membrane", membrane_gids, membrane_bound)
    ];    
    ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    push!(ac_flux,ec_sol.fluxes["EX_30089"]/ec_sol.fluxes["EX_15903"])
    push!(membrane_capacity,ec_sol.gene_product_capacity["membrane"]/membrane_bound*100)
end

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

f = Figure(; size=(10cm, 6cm))#, backgroundcolor=:transparent)

ax = Axis(
    f[1,1];
    backgroundcolor=:transparent,
    xlabel = "Membrane bound (mg/gDW)",
    ylabel = "Acetate molecules excreted per glucose consumes",
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    ygridvisible=false,
    xgridvisible=false,
)
ax2 = Axis(
    f[1, 1], 
    yaxisposition = :right,
    ylabel = "Membrane capacity used (%)",
    ylabelsize=6pt,
    yticklabelsize=5pt,
    ygridvisible=false,
    xgridvisible=false,
)
hidespines!(ax2)
hidexdecorations!(ax2)
lines!(
    ax,
    vols,
    abs.(ac_flux),
    label = "Acetate exchange"
)
lines!(
    ax2,
    vols,
    membrane_capacity,
    label = "Membrane capacity",
    color = :red
)
axislegend(
    ax,
    position=:rt,
    labelsize = 5pt,
)
f

save("plots/acetate_exchange.png", f, px_per_unit = 1200/inch)
