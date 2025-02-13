function oxphos!(model)

    df = DataFrame(CSV.File("data/model/oxphos_reactions.csv"))
    ms = RheaReactions.RheaMetabolite[]

    for row in eachrow(unique(df))
        rid = parse(Int64, split(row.RHEA_ID, ':')[2])
        haskey(model.reactions, string(rid)) && continue

        rxn = get_reaction(rid)

        coeff_mets = get_reaction_metabolites(rid)

        stoichiometry = Dict(
            string(v.accession) => s
            for (s, v) in coeff_mets
        )

        append!(ms, last.(coeff_mets))

        ecs = isnothing(rxn.ec) ? [row.EC] : [rsplit(x, '/'; limit=2)[2] for x in rxn.ec]
        name = rxn.name

        model.reactions[string(rid)] = CM.Reaction(;
            name=name,
            lower_bound=-1000.0,
            upper_bound=1000.0,
            stoichiometry=stoichiometry,
            annotations=Dict(
                "REACTION" => [rxn.equation],
                "EC" => ecs,
                "KEGG" => [string(row.KEGG_ID)]
            ),
            notes=Dict(
                "reason" => ["oxphos"]
            ),
        )
    end

    # add metabolites 
    for m in ms
        haskey(model.metabolites, m.accession) && continue
        model.metabolites[m.accession] = CM.Metabolite(
            m.name,
            nothing,
            parse_formula(m.formula),
            m.charge,
            0.0,
            Dict{String,Vector{String}}(),
            Dict{String,Vector{String}}(),
        )
    end

    model
end