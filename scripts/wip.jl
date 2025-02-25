using StaphylococcusAureus
using COBREXA, JSONFBCModels
import AbstractFBCModels.CanonicalModel as CM
using HiGHS, JSON

model = build_model()

save_model(convert(JSONFBCModels.JSONFBCModel, model), "data/model.json")

# make model with gene ids as reaction names
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")

id_tag("lcl|AM990992.1_prot_CAQ50558.1_2118")
####################################

fba_sol = parsimonious_flux_balance_analysis(model; optimizer=HiGHS.Optimizer)

open("data/fluxes.json", "w") do io
    JSON.print(io, Dict(string(x) => y for (x, y) in fba_sol.fluxes))
end

ex_fluxes = Dict(
    ("CHEBI:$(split(string(x),"EX_")[2])", model.metabolites["CHEBI:$(split(string(x),"EX_")[2])"].name,x) => y
    for (x, y) in fba_sol.fluxes if startswith(string(x), "EX") && y<-1e-5
)

["CHEBI:$(split(string(x),"EX_")[2]),$(model.metabolites["CHEBI:$(split(string(x),"EX_")[2])"].name)" for (x, y) in fba_sol.fluxes if startswith(string(x), "EX") && y<-1e-5]

Dict(x=>y for (x,y) in fba_sol.fluxes if abs(y)>1e-5)

open("data/fluxes.json", "w") do io
    JSON.print(io, Dict(string(x) => y for (x, y) in fba_sol.fluxes))
end


### make sinks 
for (m, met) in model.metabolites
    haskey(model.reactions, "EX_$(split(m,':')[2])") && continue
    model.reactions["EX_$(split(m,':')[2])"] = CM.Reaction(; stoichiometry=Dict(m => 1), notes=Dict("reason" => ["added as sink"]), upper_bound=0.0)
end


fba_sol = parsimonious_flux_balance_analysis(model; optimizer=HiGHS.Optimizer)





### delete unneeded sinks 

for (m, met) in model.metabolites
    if abs(fba_sol.fluxes["EX_$(split(m,':')[2])"]) < 1e-5
        delete!(model.reactions, "EX_$(split(m,':')[2])")
    end
end


fba_sol = parsimonious_flux_balance_analysis(model; optimizer=HiGHS.Optimizer)




#### get reactions that produce CoA CHEBI:57287

[(r,rxn.annotations["EC"],rxn.annotations["KEGG"]) for (r,rxn) in model.reactions if fba_sol.fluxes[r]>1e-5 && haskey(rxn.stoichiometry,"CHEBI:57287") && rxn.stoichiometry["CHEBI:57287"]>0]

[(r,rxn.annotations["EC"],rxn.annotations["KEGG"]) for (r,rxn) in model.reactions if fba_sol.fluxes[r]<-1e-5 && haskey(rxn.stoichiometry,"CHEBI:57287") && rxn.stoichiometry["CHEBI:57287"]<0]
