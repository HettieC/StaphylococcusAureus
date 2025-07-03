function build_model()

    df = DataFrame(CSV.File("data/model/reactions/metabolic_reactions.csv"))

    heteros = @rsubset(df, !iszero(:Isozyme))
    gheteros = groupby(heteros, [:RHEA_ID, :Isozyme])

    homos = @rsubset(df, iszero(:Isozyme))
    ghomos = groupby(homos, [:RHEA_ID, :Protein])

    # Build model

    model = CM.Model()

    extend_model!(model, ghomos)
    extend_model!(model, gheteros)
    # add new reactions 
    df = DataFrame(CSV.File("data/model/reactions/new_reactions.csv"))
    heteros = @rsubset(df, !iszero(:Isozyme))
    gheteros = groupby(heteros, [:RHEA_ID, :Isozyme])
    homos = @rsubset(df, iszero(:Isozyme))
    ghomos = groupby(homos, [:RHEA_ID, :Protein])
    extend_model!(model, ghomos)
    extend_model!(model, gheteros)
    gapfill!(model)
    add_sources!(model)
    add_sinks!(model)
    reaction_isozymes, kcat_dict = get_reaction_isozymes()
    add_periplasm_transporters!(model)
    add_membrane_transporters!(model,reaction_isozymes,kcat_dict)
    
    curate!(model)
    change_bounds!(model)
    add_electron_transport_chain!(model)
    add_special_isozymes!(reaction_isozymes,kcat_dict,model)
    add_fake_isozymes!(model,reaction_isozymes)
    add_genes!(model)
    add_names_pathways!(model)
    scale_biomass!(model)
    return model, reaction_isozymes
end
export build_model
