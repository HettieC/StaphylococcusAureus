using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using HiGHS

model, reaction_isozymes = build_model()

blocked_rxns = String[]
for (r,rxn) in model.reactions
    # if rxn.upper_bound == 0 
    #     rxn.upper_bound = -0.01
    #     sol = flux_balance_analysis(model; optimizer=HiGHS.Optimizer)
    #     if isnothing(sol) || sol.objective < 1e-5 
    #         push!(blocked_rxns, r)
    #         rxn.upper_bound = 0
    #     end
    # if rxn.lower_bound == 0 
    #     rxn.lower_bound = 0.01
    #     sol = flux_balance_analysis(model; optimizer=HiGHS.Optimizer)
    #     if isnothing(sol) || sol.objective < 1e-5 
    #         push!(blocked_rxns, r)
    #         rxn.lower_bound = 0
    #     end
    if rxn.lower_bound == -1000 && rxn.upper_bound == 1000
        rxn.lower_bound = 0.01
        sol = flux_balance_analysis(model; optimizer=HiGHS.Optimizer)
        if isnothing(sol) || sol.objective < 1e-5 
            push!(blocked_rxns, r)
            rxn.lower_bound = -1000
        end
        rxn.lower_bound = -1000
        rxn.upper_bound = -0.01 
        sol = flux_balance_analysis(model; optimizer=HiGHS.Optimizer)
        if isnothing(sol) || sol.objective < 1e-5 
            push!(blocked_rxns, r)
            rxn.upper_bound = 1000
        end
        rxn.upper_bound = 1000
    end
end
unique(blocked_rxns)

