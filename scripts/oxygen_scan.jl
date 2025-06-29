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
model.reactions["EX_15379"].upper_bound = 1000
model.reactions["EX_15379"].lower_bound = 0
model.reactions["EX_15378"].upper_bound = 0

sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
open("data/fluxes.json","w") do io 
    JSON.print(io,sol.fluxes)
end
sol.fluxes["EX_15379"]

# scan growth rate across limited oxygen with fba model 
growth = Float64[]
o2_iter = 0:1:40
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

#o2 producing 
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:15379") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:15379"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:15379"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)

# CHEBI:28938 nh4+ producing 
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:28938") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:28938"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:28938"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)

# oxygen scan with EC model 
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)
model.reactions["EX_16236"].lower_bound = 0 #block ethanol exchange
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_16651"].lower_bound = 0 #block (s)-lactate exchange
model.reactions["EX_16004"].lower_bound = 0 #block (r)-lactate exchange
model.reactions["EX_15903"].upper_bound = 1000 #unlimit glucose
model.reactions["EX_15379"].lower_bound = 0 #unlimit o2
model.reactions["EX_15379"].upper_bound = 1000


capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 180.0),
    ("membrane", membrane_gids, 60.0)
];

ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
) 
ec_sol.fluxes["EX_15379"]

ec_growth = Float64[]
for o2_uptake in o2_iter
    model.reactions["EX_15379"].upper_bound = o2_uptake
    model.reactions["EX_15379"].lower_bound = o2_uptake-0.1
    ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )    
    push!(ec_growth,ec_sol.objective)
end
ec_growth

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

f = Figure(; size=(10cm, 6cm))#, backgroundcolor=:transparent)

ax = Axis(
    f[1,1];
    backgroundcolor=:transparent,
    ylabel = "Growth rate (1/h)",
    xlabel = "Oxygen uptake rate (mMol/h)",
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    ygridvisible=false,
    xgridvisible=false,
)
lines!(
    ax,
    o2_iter,
    ec_growth,
    label = "ecFBA",
)
lines!(
    ax,
    o2_iter,
    growth,
    label = "FBA",
    color=:red,
    linestyle=:dash
)
axislegend(
    ax,
    position=:rb,
    labelsize = 5pt,
)
f
save("data/plots/o2_scan.png", f, px_per_unit = 1200/inch)

C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

