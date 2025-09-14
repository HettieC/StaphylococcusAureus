using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using HiGHS, JSON, CSV
using JSONFBCModels, DataFrames
using Latexify, CairoMakie

model, reaction_isozymes = build_model();
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)

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
    "15361" => "pyruvate",
    "16651" => "L-lactate"
)

model.reactions["EX_15903"].upper_bound = 0
model.reactions["EX_47013"].upper_bound = 0
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange

capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 200.0),
    ("membrane", membrane_gids, 130.0)
];

biomass_carbon = 0 
for (x,y) in model.reactions["biomass"].stoichiometry 
    biomass_carbon -= haskey(A.metabolite_formula(model,x),"C") ? A.metabolite_formula(model,x)["C"]*y : 0 
end
biomass_carbon

fba_growth = Dict{String,Vector{Any}}()
for (x,y) in carbon_sources 
    "CHEBI:$x" ∉ A.metabolites(model) && continue

    if !haskey(model.reactions,"EX_$x")
        println("rxn added: $y\n")
        model.reactions["EX_$x"] = CM.Reaction(;
            name = "$y exchange into cytosol",
            stoichiometry = Dict("CHEBI:$x" => 1),
            lower_bound = 0.0
        )
    end
    model.reactions["EX_$x"].upper_bound = 1#/(chebi_mass["CHEBI:$x"]/1000)
    sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    fba_growth[y] = isnothing(sol) ? ["CHEBI:$x",0] : ["CHEBI:$x",sol.objective]
    fba_growth[y] == 0 && continue

    model.reactions["EX_$x"].upper_bound = 0

end
fba_growth


# scale growth to gDW/g substrate
# get molar mass of substrates 
chebi_mass = Dict{String,Float64}()
open("data/databases/chebi/chebi_core.obo","r") do io 
    i = 0
    chebi = "" 
    mass = 0
    for ln in eachline(io)
        i += 1
        i < 20 && continue
        if startswith(ln,"id: ")
            chebi = string(split(ln,"id: ")[2])
        elseif startswith(ln,"property_value: http://purl.obolibrary.org/obo/chebi/mass ")
            mass = parse(Float64,split(ln)[3][2:end-1])    
        elseif ln == "[Term]"
            chebi_mass[chebi] = mass
        end
    end
end

c_per_c = Dict(
    x => biomass_carbon * y[2] / A.metabolite_formula(model,y[1])["C"]
    for (x,y) in fba_growth
)


df = DataFrame(Carbon_source=String[],Growth=Float64[])
for (x,y) in c_per_c
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
