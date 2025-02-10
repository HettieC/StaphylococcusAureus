using CSV, DataFrames, DataFramesMeta
using COBREXA, AbstractFBCModels
import AbstractFBCModels.CanonicalModel as CM
import COBREXA as X
using DocStringExtensions
using RheaReactions
using JSONFBCModels

include("src/utils.jl")
include("src/reconstruct.jl")

model = build_model()

save_model(convert(JSONFBCModels.JSONFBCModel,model),"data/model.json")

df = DataFrame(CSV.File("data/model/gene_ids.csv"))
g_model = deepcopy(model)
for row in eachrow(df)
    g_model.reactions[string(row.RHEA_ID)].name = row.Gene_ID 
end

save_model(convert(JSONFBCModels.JSONFBCModel,g_model),"escher_model.json")
