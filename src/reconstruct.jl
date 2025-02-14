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
    change_bounds!(model)

    model
end
export build_model
