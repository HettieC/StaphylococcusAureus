using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
using COBREXA, DataFrames, CSV
using HiGHS, JSON, DataFramesMeta
using JSONFBCModels, RheaReactions
import ConstraintTrees as C

# add this to transporters.csv: Permease,glucose,CHEBI:15903,SAPIG2309,1
model, reaction_isozymes = build_model()

biocyc = DataFrame(CSV.File(joinpath("data", "databases", "rhea", "biocyc_rxns.csv")))
bidirectional = string.(JSON.parsefile("data/model/bidirectional.json"))
@select!(biocyc, :rheaDir, :metacyc)
for rid in A.reactions(model)
    rid ∈ bidirectional && continue
    isnothing(tryparse(Int,rid)) && continue
    qrt = RheaReactions.get_reaction_quartet(parse(Int, rid))
    df = @subset(biocyc, in.(:rheaDir, Ref(qrt)))
    isempty(df) && continue
    lb, ub = rhea_rxn_dir(df[1, 1], qrt)
    model.reactions[rid].lower_bound = lb
    model.reactions[rid].upper_bound = ub
end






escher_model = change_reaction_names(model)
df = DataFrame(CSV.File("data/model/new_reactions.csv"))
for row in eachrow(df) 
    if haskey(escher_model.reactions,string(split(row.RHEA_ID,':')[2]))
        escher_model.reactions[string(split(row.RHEA_ID,':')[2])].name *= " new"
    end
end
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")

model.reactions["EX_16236"].lower_bound = 0 #block ethanol exchange
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_15903"].upper_bound = 10 #glucose bound

sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)

open("data/fluxes.json","w") do io 
    JSON.print(io,sol.fluxes)
end
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)
C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 500    
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)



####### close exchanges and see if atp can be produced 
model, reaction_isozymes = build_model()
for (r,rxn) in model.reactions 
    if startswith(r,"EX")
        rxn.upper_bound = 0 
        rxn.lower_bound = 0 
    end
end
model.reactions["biomass"].objective_coefficient = 1
model.reactions["ATPS"].objective_coefficient = 1
sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)

open("data/fluxes.json","w") do io 
    JSON.print(io,sol.fluxes)
end
