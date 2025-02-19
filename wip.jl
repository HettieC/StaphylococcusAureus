using StaphylococcusAureus
using COBREXA, JSONFBCModels
import AbstractFBCModels.CanonicalModel as CM
using HiGHS, JSON

model = build_model()

save_model(convert(JSONFBCModels.JSONFBCModel, model), "data/model.json")

# make model with gene ids as reaction names
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")

id_tag("lcl|AM990992.1_prot_CAQ48875.1_435")
####################################

model = build_model()
model.reactions["biomass"] = CM.Reaction(
    ;
    name="Biomass based on Staph epidermis RP62A",
    lower_bound=0.0,
    upper_bound=1000.0,
    stoichiometry=Dict(
        "CHEBI:30616" => -1, #atp
        "CHEBI:15377" => -1, #h2o
        "CHEBI:43474" => 1, #phosphate
        "CHEBI:15378" => 1, #h+
        "CHEBI:456216" => 1, #adp
        "CHEBI:46398" => -0.1,   #UTP
        "CHEBI:37565" => -0.059,   #GTP
        "CHEBI:37563" => -0.059,   #CTP       
        "CHEBI:57692" => -0.007,    #FAD  
        "CHEBI:61404" => -0.02,    #dATP
        "CHEBI:57287" => -4.42e-5, #CoA
        "CHEBI:37568" => -0.02,    #dTTP
        "CHEBI:61429" => -0.099,   #dGTP
        "CHEBI:61481" => -0.099,   #dCTP
        #
        #"CHEBI:57783" => -2e-5, #NADPH 
        #"CHEBI:57945" => -2e-5, #NADH
        #"CHEBI:58349" => -2e-5, #NADP(+)
        #"CHEBI:57540" => -2e-5, #NAD(+)


        #"CHEBI:30807" => -1.0,    #tetradecanoate
        "CHEBI:25646" => -1.0,    #octanoate
        #"CHEBI:7896" => -1.0,     #hexadecanoate
        #"CHEBI:18262" => -1.0,    #dodecanoate
        "CHEBI:27689" => -1.0,    #decanoate
        "CHEBI:57427" => -0.282,  #L-leucine
        "CHEBI:32682" => -0.111,  #L-arginine  
        "CHEBI:57762" => -0.207,  #L-valine  
        "CHEBI:60039" => -0.116,  #L-proline
        "CHEBI:35235" => -0.019,  #L-cysteine
        "CHEBI:57305" => -0.19,  #glycine
        "CHEBI:33384" => -0.198, #L-serine         
        "CHEBI:29991" => -0.261, #L-aspartate
        "CHEBI:57972" => -0.212, #L-alanine
        "CHEBI:58359" => -0.2,   #L-glutamine
        "CHEBI:57912" => -1.0,   #L-tryptophan
        "CHEBI:29985" => -1.0,   #L-glutamate
        "CHEBI:32551" => -0.235, #L-lysine
        "CHEBI:57844" => -0.084, #L-methionine
        "CHEBI:58045" => -0.269, #L-isoleucine
        "CHEBI:57595" => -0.073, #L-histidine
        "CHEBI:57926" => -0.179, #L-threonine
        "CHEBI:58095" => -0.137, #L-phenylalanine
        "CHEBI:58315" => -0.119, #L-tyrosine
        "CHEBI:58048" => -0.1,  #L-asparagine
        "CHEBI:57305" => -0.1,  #glycine
        #"CHEBI:59789" => -0.084,
        #"CHEBI:57856" => -0.084,
    ),
    objective_coefficient=1.0,
    notes=Dict("ref" => ["Diaz Calvo, S. epidermis, Metabolites 2022"]),
)

fba_sol = flux_balance_analysis(model; optimizer=HiGHS.Optimizer)

ex_fluxes = Dict(
    ("CHEBI:$(split(string(x),"EX_")[2])", model.metabolites["CHEBI:$(split(string(x),"EX_")[2])"].name,x) => y
    for (x, y) in fba_sol.fluxes if startswith(string(x), "EX")
)

open("fluxes.json", "w") do io
    JSON.print(io, Dict(string(x) => y for (x, y) in fba_sol.fluxes))
end

open("fluxes_rhea.json", "w") do io
    JSON.print(io, Dict("RHEA:$(string(x))" => y for (x, y) in fba_sol.fluxes))
end


delete!(model.reactions,"EX_29969")

