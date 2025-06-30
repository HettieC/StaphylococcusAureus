using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using HiGHS, JSON, CSV
using JSONFBCModels, DataFrames
using Latexify

model, reaction_isozymes = build_model()

carbon_sources = Dict(
    "15903" => "beta-D-glucose",
    "47013" => "D-ribose",
    "57972" => "alanine",
    "30031" => "succinate",
    "17754" => "glycerol",
    "29985" => "glutamate",
    "47013" => "ribose",
    "30089" => "acetate",
    "58723" => "glucosamine",
    "506227" => "n-acetyl-d-glucosamine",
    "18391" => "gluconate",
    "15589" => "malate",
    "29806" => "fumarate",
    "40886" => "L-arabinose",
    "15361" => "pyruvate"
)

model.reactions["EX_15903"].upper_bound = 0
model.reactions["EX_47013"].upper_bound = 0

growth = Dict{String,Float64}()
for (x,y) in carbon_sources 
    "CHEBI:$x" ∉ A.metabolites(model) && continue
    println("\nmetab: $y")

    if !haskey(model.reactions,"EX_$x")
        println("rxn added \n")
        model.reactions["EX_$x"] = CM.Reaction(;
            name = "$y exchange in cytosol",
            stoichiometry = Dict("CHEBI:$x" => 1),
            lower_bound = 0.0
        )
    end
    model.reactions["EX_$x"].upper_bound = 10
    sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    delete!(model.reactions,"EX_$x")
    isnothing(sol) && continue 
    # C.pretty(
    #     C.ifilter_leaves(sol.fluxes) do ix, x
    #         abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    #     end; 
    #     format_label = x -> A.reaction_name(model, string(last(x))),
    # )
    growth[y] = sol.objective 
    println(sol.fluxes["EX_$x"])
end
growth

df = DataFrame(Carbon_source=String[],Growth=Float64[])
for (x,y) in growth
    push!(
            df,
            [
                x,
                round(y,sigdigits=4),
            ])
end
df
sort!(df,:Growth)
latexify(df; env = :table, booktabs = true, latex = false) |> print

#
