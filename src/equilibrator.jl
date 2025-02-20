### only needs to be done once

using PyCall
using Conda

Conda.add("equilibrator-api"; channel="conda-forge") # note: use a dash and not an underscore in "equilibrator-api"
eq = pyimport("equilibrator_api") # note: use an underscore and not a dash in "equilibrator_api"

using eQuilibrator, Unitful, CSV, DataFrames, DataFramesMeta, RheaReactions, Measurements

df = DataFrame(CSV.File("data/model/metabolic_reactions.csv"))

eq = eQuilibrator.Equilibrator(); # slow, only do if necessary

ds = Union{Missing, Float64}[]
pdgs = Union{Missing, Float64}[]

for r in eachrow(df)
    println(r)
    
    rid = parse(Int64,split(r.RHEA_ID,':')[2])
    coeff_mets = get_reaction_metabolites(rid)

    stoichiometry = Dict(
        string(v.accession) => s
        for (s, v) in coeff_mets
    )

    substrates = [string(abs(Int(v))) * " chebi:" * k for (k, v) in stoichiometry if v < 0]
    products = [string(abs(Int(v))) * " chebi:" * k for (k, v) in stoichiometry if v > 0]
    rxn_string = join(substrates, " + ") * " = " * join(products, " + ")
    
    try 
        _d = ln_reversibility_index(eq, rxn_string; skip_unbalanced = true)
        _pdg =  physiological_dg_prime(eq, rxn_string; skip_unbalanced = true)
        d = isnothing(_d) ? missing : Measurements.value(_d)
        pdg = isnothing(_pdg) ? missing : Measurements.value(ustrip(u"kJ/mol", _pdg))
    catch
        d = missing
        pdg = missing
    end

    push!(ds, d)
    push!(pdgs, pdg)
end

df.DeltaG = pdgs
df.RevIndex = ds
df

CSV.write("data/model/metabolic_reactions.csv",df)
