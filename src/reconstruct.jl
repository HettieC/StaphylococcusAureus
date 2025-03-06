function build_model()

    df = DataFrame(CSV.File("data/model/metabolic_reactions.csv"))

    heteros = @rsubset(df, !iszero(:Subunit))
    gheteros = groupby(heteros, [:RHEA_ID, :Subunit])

    homos = @rsubset(df, iszero(:Subunit))
    ghomos = groupby(homos, [:RHEA_ID, :Protein])

    # Build model

    model = CM.Model()

    extend_model!(model, ghomos)
    extend_model!(model, gheteros)
    gapfill!(model)
    add_sources!(model)
    add_sinks!(model)
    add_oxphos!(model)
    model = curate!(model)
    change_bounds!(model)
    add_electron_transport_chain!(model)

    # remove general quinone/quinol/ubiquinone/ubiquinone reactions 
    for (r,rxn) in model.reactions 
        if haskey(rxn.stoichiometry,"CHEBI:132124") || haskey(rxn.stoichiometry,"CHEBI:24646") || haskey(rxn.stoichiometry,"POLYMER:9566") || haskey(rxn.stoichiometry,"POLYMER:9565") || haskey(rxn.stoichiometry,"POLYMER:9563")
            delete!(model.reactions,r)
        end
    end

    #

    return model
end
export build_model
