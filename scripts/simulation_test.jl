using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
using COBREXA, JSONFBCModels
using HiGHS, JSON
using ConstraintTrees
import ConstraintTrees as C

model = build_model()

save_model(convert(JSONFBCModels.JSONFBCModel, model), "data/model.json")
# make model with gene ids as reaction names
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")

model.reactions["EX_15903"].upper_bound = 10 #glucose
sol = flux_balance_analysis(model, optimizer=HiGHS.Optimizer)

#######
### need some arginine symport from periplasm!

#######

sol = parsimonious_flux_balance_analysis(model, optimizer=HiGHS.Optimizer)

open("data/fluxes.json","w") do io 
    JSON.print(io,sol.fluxes)
end

open("data/big_fluxes.json","w") do io 
    JSON.print(io,Dict(x=>y for (x,y) in sol.fluxes if abs(y)>50))
end

Dict((x,model.reactions[String(x)].name)=>y for (x,y) in sol.fluxes if abs(y)>15)


# add sinks
for (m, met) in model.metabolites
    haskey(model.reactions, "EX_$(split(m,':')[2])") && continue
    model.reactions["EX_$(split(m,':')[2])"] = CM.Reaction(; stoichiometry=Dict(m => 1), notes=Dict("reason" => ["added as sink"]), upper_bound=0.0)
end


sol = parsimonious_flux_balance_analysis(model; optimizer=HiGHS.Optimizer)

### delete unneeded sinks 

for (m, met) in model.metabolites
    if abs(fba_sol.fluxes["EX_$(split(m,':')[2])"]) < 1e-5
        delete!(model.reactions, "EX_$(split(m,':')[2])")
    end
end


sol = parsimonious_flux_balance_analysis(model; optimizer=HiGHS.Optimizer)


open("data/fluxes.json", "w") do io
    JSON.print(io, Dict(string(x) => y for (x, y) in sol.fluxes))
end

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && isnothing(tryparse(Int,string(last(ix))))    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

C.pretty(
    C.ifilter_leaves(loopless_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

# atp producing reactions
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:30616") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)



# proton producing reactions
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:15378") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:15378"] > 0 && x > 1) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:15378"] < 0 && x < -1))
    end; 
)
