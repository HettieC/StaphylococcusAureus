using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
using COBREXA, DataFrames, CSV
using HiGHS, JSON, DataFramesMeta
using JSONFBCModels, RheaReactions
import ConstraintTrees as C

# add this to transporters.csv: Permease,glucose,CHEBI:15903,SAPIG2309,1
model, reaction_isozymes = build_model()


sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)


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


sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)

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
model.reactions["biomass"].objective_coefficient = 0
model.reactions["ATPS"].objective_coefficient = 0

model.reactions["ATPM"].objective_coefficient = 1
sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)

open("data/fluxes.json","w") do io 
    JSON.print(io,sol.fluxes)
end


using RheaReactions, DataFramesMeta
    #change directions to match what is found in biocyc 
biocyc = DataFrame(CSV.File(joinpath("data", "databases", "rhea", "biocyc_rxns.csv")))
manual_directions = DataFrame(CSV.File("data/model/unidirectional_reactions.csv"))
#bidirectional = string.(JSON.parsefile("data/model/bidirectional.json"))
@select!(biocyc, :rheaDir, :metacyc)
bidirectional = String[]
for rid in A.reactions(model)
    #rid ∈ bidirectional && continue
    isnothing(tryparse(Int,rid)) && continue
    rid ∈ string.(manual_directions.RHEA_ID) && continue
    qrt = RheaReactions.get_reaction_quartet(parse(Int, rid))
    df = @subset(biocyc, in.(:rheaDir, Ref(qrt)))
    isempty(df) && continue
    lb, ub = rhea_rxn_dir(df[1, 1], qrt)
    model.reactions[rid].lower_bound = lb
    model.reactions[rid].upper_bound = ub
    sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    if isnothing(sol) || sol.objective < 1e-3 
        push!(bidirectional,rid)
        model.reactions[rid].lower_bound = -1000
        model.reactions[rid].upper_bound = 1000
    end
end

open("data/model/bidirectional.json","w") do io 
    JSON.print(io,bidirectional)
end

[r for r in manual_directions.RHEA_ID if haskey(model.reactions,string(r)) && model.reactions[string(r)].lower_bound == -1000 && model.reactions[string(r)].upper_bound == 1000]




directions = Dict(
    rid => [rxn.lower_bound,rxn.upper_bound] for (rid,rxn) in model.reactions
)
open("data/model/directions2.json","w") do io 
    JSON.print(io,directions)
end

dir = Dict(x=>y for (x,y) in JSON.parsefile("data/model/directions.json"))

Dict(x=>(y,dir[x]) for (x,y) in directions if y != dir[x])
