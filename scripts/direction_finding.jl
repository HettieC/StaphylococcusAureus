using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using DataFrames, CSV
using HiGHS, DataFramesMeta
using JSON, JSONFBCModels, RheaReactions

# add this to transporters.csv: Permease,glucose,CHEBI:15903,SAPIG2309,1

model, reaction_isozymes = build_model()

# check respiration 
exch_ids = ["EX_15377","EX_15379","EX_15903","EX_16526"]
for rid in [r for r in A.reactions(model) if startswith(r,"EX")]
    if rid ∉ exch_ids 
        model.reactions[rid].lower_bound = model.reactions[rid].upper_bound = 0 
    end
end
model.reactions["EX_15903"].upper_bound = 1
model.reactions["biomass"].stoichiometry = Dict(
            "CHEBI:15903" => -1, #glucose
            "CHEBI:15377" => 6, #h2o
            "CHEBI:15379" => -6, #o2
            "CHEBI:16526" => 6, #co2
        )
ec_sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

model.reactions["biomass"].stoichiometry = Dict(
            "CHEBI:30616" => -50, #atp
            "CHEBI:15377" => -50, #h2o
            "CHEBI:43474" => 50, #phosphate
            "CHEBI:15378" => 50, #h+
            "CHEBI:456216" => 50, #adp
        )
ec_sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
ec_sol.fluxes["EX_30089"] #acetate 
# make model with gene ids as reaction names
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")

open("data/fluxes.json","w") do io 
    JSON.print(io,ec_sol.fluxes)
end

# atp producing reactions
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:30616") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)



model, reaction_isozymes = build_model()
df = DataFrame(CSV.File("data/model/unidirectional_reactions.csv"))
model.reactions["EX_15903"].upper_bound = 10 #glucose
model.reactions["EX_47013"].upper_bound = 0 #ribose
rxns = []
for ln in eachrow(df) 
    !haskey(model.reactions,string(ln.RHEA_ID)) && continue
    model.reactions[string(ln.RHEA_ID)].lower_bound = ln.LOWER_BOUND + 0.0 
    model.reactions[string(ln.RHEA_ID)].upper_bound = ln.UPPER_BOUND + 0.0
    ec_sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    if isnothing(ec_sol) || abs(ec_sol.objective) < 1e-5 || -10000 < ec_sol.fluxes["EX_30089"] < -1e-5
        push!(rxns,ln.RHEA_ID)
        model.reactions[string(ln.RHEA_ID)].lower_bound = -1000 
        model.reactions[string(ln.RHEA_ID)].upper_bound = 1000
    end
end
ec_sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
ec_sol.fluxes["EX_30089"] #acetate 
# make model with gene ids as reaction names
escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")

open("data/fluxes.json","w") do io 
    JSON.print(io,ec_sol.fluxes)
end
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 100    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)


C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-6 && haskey(A.reaction_stoichiometry(model, string(first(ix))),"CHEBI:30616") && ( (model.reactions[string(first(ix))].stoichiometry["CHEBI:30616"]>0 && x>1e-5) || (model.reactions[string(first(ix))].stoichiometry["CHEBI:30616"]<0 && x<1e-5) )
    end; 
    #format_label = x -> A.reaction_name(model, string(last(x))),
)

rhea_rxn_dir(rxn, qrt) = begin
    idx = first(indexin([rxn], qrt))
    isnothing(idx) && error("Reaction not found...")
    idx == 1 && return (-1000, 1000)
    idx == 2 && return (0, 1000)
    idx == 3 && return (-1000, 0)
    idx == 4 && return (-1000, 1000)
end
# change directions to match what is found in biocyc - manual thermodynamics leaves much to be desired
biocyc = DataFrame(CSV.File(joinpath("data", "databases", "rhea", "biocyc_rxns.csv")))
@select!(biocyc, :rheaDir, :metacyc)
directions = String[]
for rid in A.reactions(model)
    isnothing(tryparse(Int,rid)) && continue
    qrt = RheaReactions.get_reaction_quartet(parse(Int, rid))
    df = @subset(biocyc, in.(:rheaDir, Ref(qrt)))
    isempty(df) && continue
    lb, ub = rhea_rxn_dir(df[1, 1], qrt)
    model.reactions[rid].lower_bound = lb
    model.reactions[rid].upper_bound = ub
    if lb != -1000 || ub != 1000
        push!(directions, rid)
    end
end
