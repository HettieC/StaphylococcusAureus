using CSV, DataFrames, DataFramesMeta
using COBREXA, AbstractFBCModels
import AbstractFBCModels.CanonicalModel as CM
import COBREXA as X
using DocStringExtensions
using RheaReactions
using JSONFBCModels, JSON
using HiGHS

include("src/utils.jl")
include("src/reconstruct.jl")

model = build_model()

save_model(convert(JSONFBCModels.JSONFBCModel,model),"data/model.json")

# make model with gene ids as reaction names
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel,escher_model),"escher_model.json")


####################################


model.reactions["biomass"] = CM.Reaction(;
    name = "biomass",
    lower_bound = 0.0,
    upper_bound = 1000.0,
    stoichiometry = Dict(
        "CHEBI:30616" => -1, #atp
        "CHEBI:15377" => -1, #h2o
        "CHEBI:43474" => 1, #phosphate
        "CHEBI:15378" => 1, #h+
        "CHEBI:456216" => 1, #adp
    ),
    objective_coefficient = 1.0
)

fba_sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)

open("fluxes.json","w") do io 
    JSON.print(io,Dict("RHEA:$(string(x))" => y for (x,y) in fba_sol.fluxes))
end
