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

save_model(convert(JSONFBCModels.JSONFBCModel, model), "data/model.json")

# make model with gene ids as reaction names
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")


####################################


model.reactions["biomass"] = CM.Reaction(;
    name="biomass",
    lower_bound=0.0,
    upper_bound=1000.0,
    stoichiometry=Dict(
        "CHEBI:30616" => -1, #atp
        "CHEBI:15377" => -1, #h2o
        "CHEBI:43474" => 1, #phosphate
        "CHEBI:15378" => 1, #h+
        "CHEBI:456216" => 1, #adp

        #"CHEBI:46398" => -0.061,   #UTP
        #"CHEBI:37565" => -0.059,  #GTP
        #"CHEBI:37563" => -0.059,   #CTP       

        #"CHEBI:57540" => -0.00158, #NAD(+)
        #"CHEBI:57783" => -2.95e-5, #NADPH       
        #"CHEBI:57945" => -3.65e-5, #NADH
        #"CHEBI:57692" => -74.1,    #FAD  
        #"CHEBI:57287" => -4.42e-5, #CoA
        #"CHEBI:58349" => -9.62e-5,#NADP(+)   

        #"CHEBI:61429" => -0.099,  #dGTP
        #"CHEBI:61481" => -0.099,  #dCTP
        #"CHEBI:37568" => -0.02,    #dTTP
        #"CHEBI:61404" => -0.02,    #dATP

        #"CHEBI:30807" => -1.0,     #tetradecanoate
        #"CHEBI:25646" => -1.0,     #octanoate
        #"CHEBI:7896" => -1.0,     #hexadecanoate
        #"CHEBI:18262" => -1.0,    #dodecanoate
        #"CHEBI:27689" => -1.0,    #decanoate

        #"CHEBI:57305" => -0.19,    #glycine
        #"CHEBI:33384" => -0.198,   #L-serine
        #"CHEBI:57427" => -0.282,   #L-leucine
        #"CHEBI:57972" => -0.212,   #L-alanine
        #"CHEBI:32682" => -0.111,  #L-arginine          
        #"CHEBI:29991" => -0.261,   #L-aspartate
        #"CHEBI:57762" => -0.207,   #L-valine     
        #"CHEBI:58359" => -0.2,     #L-glutamine
        #"CHEBI:57912" => -1.0,    #L-tryptophan
        #"CHEBI:29985" => -1.0,     #L-glutamate
        #"CHEBI:32551" => -0.235,   #L-lysine
        #"CHEBI:57844" => -0.084,   #L-methionine
        #"CHEBI:60039" => -0.116,   #L-proline
        #"CHEBI:58045" => -0.269,   #L-isoleucine
        #"CHEBI:57595" => -0.073,    #L-histidine
        #"CHEBI:57926" => -0.179,    #L-threonine
        #"CHEBI:35235" => -0.019,   #L-cysteine
        #"CHEBI:58095" => -0.137,   #L-phenylalanine
        #"CHEBI:58315" => -0.119,   #L-tyrosine
    ),
    objective_coefficient=1.0
)

fba_sol = flux_balance_analysis(model; optimizer=HiGHS.Optimizer)

open("fluxes.json", "w") do io
    JSON.print(io, Dict("RHEA:$(string(x))" => y for (x, y) in fba_sol.fluxes))
end
