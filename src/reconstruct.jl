function build_model()

    df = DataFrame(CSV.File("data/model/metabolic_reactions.csv"))

    heteros = @rsubset(df, !iszero(:Isozyme))
    gheteros = groupby(heteros, [:RHEA_ID, :Isozyme])

    homos = @rsubset(df, iszero(:Isozyme))
    ghomos = groupby(homos, [:RHEA_ID, :Protein])

    # Build model

    model = CM.Model()

    extend_model!(model, ghomos)
    extend_model!(model, gheteros)
    gapfill!(model)
    add_sources!(model)
    add_sinks!(model)
    reaction_isozymes, kcat_dict = get_reaction_isozymes()
    add_periplasm_transporters!(model)
    add_membrane_transporters!(model,reaction_isozymes,kcat_dict)
    add_oxphos!(model)

    model = curate!(model)
    change_bounds!(model)
    add_electron_transport_chain!(model)
    add_special_isozymes!(reaction_isozymes,kcat_dict,model)

    # remove general quinone/quinol/ubiquinone/ubiquinone reactions 
    # for (r,rxn) in model.reactions 
    #     if haskey(rxn.stoichiometry,"CHEBI:132124") || haskey(rxn.stoichiometry,"CHEBI:24646") || haskey(rxn.stoichiometry,"POLYMER:9566") || haskey(rxn.stoichiometry,"POLYMER:9565") || haskey(rxn.stoichiometry,"POLYMER:9563")
    #         delete!(model.reactions,r)
    #     end
    # end

    #

    return model, reaction_isozymes
end
export build_model
