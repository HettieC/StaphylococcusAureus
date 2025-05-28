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
