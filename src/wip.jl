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
molar_masses = Dict{String,Float64}() # g/mol
begin
    molar_masses["CHEBI:61404"] = 487.1499 # dATP
    molar_masses["CHEBI:61429"] = 503.1493 # dGTP
    molar_masses["CHEBI:61481"] = 463.1252 # dCTP
    molar_masses["CHEBI:37568"] = 478.1365 # dTTP

    molar_masses["CHEBI:30616"] = 503.14946 # ATP
    molar_masses["CHEBI:37565"] = 519.14886 # GTP
    molar_masses["CHEBI:37563"] = 479.12468 # CTP
    molar_masses["CHEBI:46398"] = 480.1094 # UTP

    molar_masses["CHEBI:32551"] = 147.19558 # # lysine
    molar_masses["CHEBI:58045"] = 131.1729 # isoleucine
    molar_masses["CHEBI:57427"] = 131.1729 # leucine
    molar_masses["CHEBI:57844"] = 149.2124 # methionine
    molar_masses["CHEBI:58095"] = 165.1891 # phenylalanine
    molar_masses["CHEBI:57926"] = 119.1197  # threonine
    molar_masses["CHEBI:57912"] = 204.2262 # tryptophan
    molar_masses["CHEBI:57762"] = 117.1469 # valine
    molar_masses["CHEBI:32682"] = 175.20906 # arginine
    molar_masses["CHEBI:57595"] = 155.1552 # histidine
    molar_masses["CHEBI:57972"] = 89.0935 # alanine
    molar_masses["CHEBI:58048"] = 132.1184 # asparagine
    molar_masses["CHEBI:29991"] = 132.09478 # aspartate
    molar_masses["CHEBI:35235"] = 121.159 # cysteine
    molar_masses["CHEBI:29985"] = 147.1299 # glutamate
    molar_masses["CHEBI:58359"] = 146.1451 # glutamine
    molar_masses["CHEBI:57305"] = 75.0669 # glycine
    molar_masses["CHEBI:60039"] = 115.131 # proline
    molar_masses["CHEBI:33384"] = 105.093 # serine
    molar_masses["CHEBI:58315"] = 181.1894 # tyrosine

    #molar_masses["glycogen"] = 162.1406 # C6H10O5

    molar_masses["CHEBI:30807"] = 227.364 # tetradecanoic acid
    molar_masses["CHEBI:7896"] = 255.4161 # hexadecanoic acid
    molar_masses["CHEBI:25629"] = 283.47 # octadecanoic acid

    #molar_masses["peptidoglycan"] = 1916.20990

    #molar_masses["kdo_lps"] = 2232.67080

    # soluble pool
    molar_masses["CHEBI:60530"] = 836.838 # Fe(II)-heme o
    molar_masses["CHEBI:57692"] = 782.5259 # FAD
    molar_masses["CHEBI:57705"] = 605.3378 # UDP-N-acetyl-alpha-D-glucosamine
    molar_masses["CHEBI:57540"] = 662.4172 # NAD(+)
    molar_masses["CHEBI:58885"] = 564.2859 # UDP-alpha-D-glucose
    molar_masses["CHEBI:57287"] = 763.502 # CoA
    molar_masses["CHEBI:57925"] = 306.31 # glutathione
    molar_masses["CHEBI:57945"] = 663.4251 # NADH
    molar_masses["CHEBI:58223"] = 401.1374 # UDP
    molar_masses["CHEBI:29985"] = 146.12136 # L-glutamate
    molar_masses["CHEBI:32966"] = 336.08392 # beta-D-fructose 1,6-bisphosphate
    molar_masses["CHEBI:30616"] = 503.14946 # ATP
    molar_masses["CHEBI:57783"] = 741.3891 # NADPH
    molar_masses["CHEBI:57986"] = 375.356 # riboflavin
    #molar_masses["CHEBI:597326"] = 245.126 # pyridoxal 5'-phosphate
    molar_masses["CHEBI:62501"] = 439.3816 # folate
    molar_masses["CHEBI:58297"] = 610.615 # glutathione disulfide
    molar_masses["CHEBI:58210"] = 453.325 # FMN
    molar_masses["CHEBI:58349"] = 740.3812 # NADP(+)
end
# convert units from g/mol to kg/mol

for (met,mass) in molar_masses
    molar_masses[met] = molar_masses[met] / 1000
end
model = build_model()
model.reactions["biomass"] = CM.Reaction(
    ;
    name="Biomass based on Staph epidermis RP62A",
    lower_bound=0.0,
    upper_bound=1000.0,
    stoichiometry = molar_masses
    ,
    objective_coefficient=1.0,
    notes=Dict("ref" => ["Diaz Calvo, S. epidermis, Metabolites 2022"]),
)

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

open("data/fluxes_rhea.json", "w") do io
    JSON.print(io, Dict("RHEA:$(string(x))" => y for (x, y) in fba_sol.fluxes))
end


delete!(model.reactions, "EX_29969")


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

