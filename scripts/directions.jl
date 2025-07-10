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

# change bounds of reactions and see if it makes nice acetate plot 
ac_flux = Dict{String,Vector{Float64}}()
biomass = Dict{String,Vector{Float64}}()
i = 0
for (r,rxn) in model.reactions
    i += 1 
    #i > 5 && break
    (rxn.lower_bound < -100 && rxn.upper_bound > 100) && continue 
    isnothing(tryparse(Int,r)) && continue 
    lb = -1000
    ub = 1000
    if rxn.lower_bound == 0 
        rxn.lower_bound = -1000 
        lb = 0
    else
        rxn.upper_bound = 1000 
        ub = 0
    end
    sol = enzyme_constrained_flux_balance_analysis(
            model;
            reaction_isozymes,
            gene_product_molar_masses,
            capacity,
            optimizer=HiGHS.Optimizer,
        )
    isnothing(sol) && continue
    abs(sol.fluxes["EX_30089"]) < 1e-4 && continue
    ac_flux[r] = [sol.fluxes["EX_30089"]]
    biomass[r] = [sol.objective]
    model.reactions["biomass"].upper_bound = copy(sol.objective)/10
    sol_new = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    push!(ac_flux[r],sol_new.fluxes["EX_30089"])
    push!(biomass[r],sol_new.objective)
    model.reactions["biomass"].upper_bound = copy(sol_new.objective)/10
    new_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    push!(ac_flux[r],new_sol.fluxes["EX_30089"])
    push!(biomass[r],new_sol.objective)
    model.reactions["biomass"].upper_bound = 1000
    model.reactions[r].lower_bound = lb 
    model.reactions[r].upper_bound = ub
end

ac_rxns = [x for (x,y) in ac_flux if abs(y[3])<1e-3 && !haskey(model.reactions[x].stoichiometry,"CHEBI:30616")]

open("maybe_bidirectional.JSON","w") do io 
    JSON.print(io,ac_rxns)
end

ac_rxns = JSON.parsefile("maybe_bidirectional.JSON")
filter!(r -> !haskey(model.reactions[r].stoichiometry,"CHEBI:30616"),ac_rxns)
inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)


#try changing the bound of a rxn and making the ac plot
i = 0
for r in ac_rxns 
    r != "10255" && continue
    i += 1
    lb = copy(model.reactions[r].lower_bound)
    ub = copy(model.reactions[r].upper_bound)

    model.reactions[r].upper_bound = 1000
    model.reactions[r].lower_bound = -1000

    ac_flux = Float64[]

    vols = 0.1:0.4:3
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
    end
    model.reactions[r].lower_bound = lb 
    model.reactions[r].upper_bound = ub

    f = Figure(; size=(10cm, 6cm))

    ax = Axis(
        f[1,1];
        backgroundcolor=:transparent,
        xlabel = "Growth rate (1/h)",
        ylabel = "Acetate exchange rate (mMol/h)",
        title = r,
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
        abs.(ac_flux),
        label = "Acetate exchange"
    )
    axislegend(
        ax,
        position=:lt,
        labelsize = 5pt,
    )
    f
    display(f)
end
