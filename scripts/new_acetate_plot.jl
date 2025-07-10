using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using CairoMakie
using HiGHS, JSON
using JSONFBCModels

model, reaction_isozymes = build_model();

model.reactions["EX_15903"].upper_bound = 10 #glucose
model.reactions["EX_47013"].upper_bound = 0 #ribose
model.reactions["EX_28938"].lower_bound = 0 #nh4+
model.reactions["EX_15378"].lower_bound = 0 #H+

fba_sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
C.pretty(
    C.ifilter_leaves(fba_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)
fba_sol.fluxes["EX_30089"]

open("data/fluxes.json","w") do io 
    JSON.print(io,fba_sol.fluxes)
end





rids = filter(x -> !startswith(x, "EX_") && x != "biomass", A.reactions(model))
unbal_rids = String[]
for rid in rids
    s = A.reaction_stoichiometry(model, rid)
    m = Dict()
    for (k, v) in s
        isnothing(A.metabolite_formula(model, k)) && continue
        for (kk, vv) in A.metabolite_formula(model, k)
            m[kk] = get(m, kk, 0) + vv * v
        end
    end
    m
    all(values(m) .== 0) || push!(unbal_rids, rid)
end
unbal_rids

rid = "10580"
s = A.reaction_stoichiometry(model, rid)
m = Dict()
for (k, v) in s
    isnothing(A.metabolite_formula(model, k)) && continue
    for (kk, vv) in A.metabolite_formula(model, k)
        m[kk] = get(m, kk, 0) + vv * v
    end
end
m



#asparagine
C.pretty(
    C.ifilter_leaves(fba_sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:58048") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:58048"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:58048"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)
Dict(r => rxn.name for (r,rxn) in model.reactions if haskey(rxn.stoichiometry,"CHEBI:58048"))

model.reactions["EX_15903"].upper_bound = 1000
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)
capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 400.0),
    ("membrane", membrane_gids, 120.0)
];
#model.reactions["EX_15378"].lower_bound = 0 #block H+ exchange
#model.reactions["EX_15740"].lower_bound = 0 #block formate exchange
#model.reactions["EX_15589"].lower_bound = 0 #block (S)-malate exchange
model.reactions["EX_16236"].lower_bound = 0 #block ethanol exchange
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
#model.reactions["EX_16651"].lower_bound = 0 #block (s)-lactate exchange
#model.reactions["EX_16004"].lower_bound = 0 #block (r)-lactate exchange
#model.reactions["EX_17544"].lower_bound = 0 #block hydrogencarbonate exchange


ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity);
ec_sol = optimized_values(
    ct;
    optimizer = HiGHS.Optimizer,
    objective = ct.objective.value,
    sense = Maximal
);
ec_sol.objective
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)


model.reactions["biomass"].upper_bound = 1
ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity);
ec_sol = optimized_values(
    ct;
    optimizer = HiGHS.Optimizer,
    objective = sum_value(ct.gene_product_amounts),
    sense = Minimal
);
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)












ac_flux = Float64[]
max_growth = copy(ec_sol.objective)
exchanges = String[]
bio_iter = 0:0.05:max_growth
for bio in bio_iter
    model.reactions["biomass"].lower_bound = bio 
    model.reactions["biomass"].upper_bound = bio+0.01
    ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity)
    
    ec_sol = optimized_values(
        ct;
        optimizer = HiGHS.Optimizer,
        objective = sum_value(ct.gene_product_amounts),
        sense = Minimal
    )
    push!(ac_flux,ec_sol.fluxes["EX_30089"])
    append!(exchanges,[string(r) for (r,v) in ec_sol.fluxes if startswith(string(r),"EX") && v<0 && string(r) ∉ exchanges])
end
ac_flux
[A.reaction_name(model,r) for r in exchanges]
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
    bio_iter,
    abs.(round.(ac_flux)),
    label = "Membrane: 120mg/gDW",
    color = :red,
    linestyle = :dash,
    linewidth=2.5

)
axislegend(
    ax,
    position=:lt,
    labelsize = 5pt,
)
display(f)


model.reactions["biomass"].lower_bound = 1 ;
model.reactions["biomass"].upper_bound = 1.01;
model.reactions["biomass"].objective_coefficient = 0
ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity);

ec_sol = optimized_values(
    ct;
    optimizer = HiGHS.Optimizer,
    objective = sum_value(ct.gene_product_capacity),
    sense = Minimal
);
ec_sol.fluxes["EX_30089"]
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)
