using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
using COBREXA, JSONFBCModels
using HiGHS, JSON
using ConstraintTrees
import ConstraintTrees as C

model = build_model()
# make model with gene ids as reaction names
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")

model.reactions["biomass"].objective_coefficient = 0.0
model.reactions["ATPM"].objective_coefficient = 1.0
model.reactions["EX_15903"].upper_bound = 1
sol = parsimonious_flux_balance_analysis(model, optimizer=HiGHS.Optimizer)

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
        abs(x) > 1e-6 && begin
            mets = [A.metabolite_name(model, k) for k in keys(A.reaction_stoichiometry(model, string(last(ix))))]
            any(in.(mets, Ref(["L-2,4-diaminobutanoate"])))
        end 
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)
