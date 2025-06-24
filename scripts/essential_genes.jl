using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using CairoMakie
using HiGHS, JSON
using JSONFBCModels, DataFrames

model, reaction_isozymes = build_model()

# save reaction names and pathways 
dic = Dict(
    r => [rxn.name, rxn.annotations["Pathway"]] for (r,rxn) in model.reactions if haskey(rxn.annotations,"Pathway")
)
open("data/model/reactions/kegg_names_pathways.json","w") do io 
    JSON.print(io,dic)
end

model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_15903"].upper_bound = 10 #limit glucose

ess_rxn = String[] 
for (r,rxn) in model.reactions 
    r == "biomass" && continue
    lb = rxn.lower_bound 
    ub = rxn.upper_bound 
    model.reactions[r].lower_bound = 0 
    model.reactions[r].upper_bound = 0 
    sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    if isnothing(sol) || sol.objective < 1e-5 
        push!(ess_rxn,r)
    end
    model.reactions[r].lower_bound = lb 
    model.reactions[r].upper_bound = ub 
end


using Latexify, KEGGReactions
df = DataFrame(ID=String[],Name=String[],Stoichiometry=String[],Pathway=String[])
for r in sort(ess_rxn)
    push!(
        df,
        [
            r,
            isnothing(model.reactions[r].name) ? "" : model.reactions[r].name ,
            !haskey(model.reactions[r].annotations,"REACTION") ? "" : model.reactions[r].annotations["REACTION"][1],
            ""
        ]
    )
end

select!(df, [:ID, :Stoichiometry])

latexify(df; env = :table, booktabs = true, latex = false) |> print

k = get_kegg_info("R00253")
