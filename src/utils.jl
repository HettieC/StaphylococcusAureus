"""
$(TYPEDSIGNATURES)

Parse the chemical formula from a RheaReaction.formula
"""
function parse_formula(x::Union{Nothing,String})
    if isnothing(x)
        return nothing
    elseif occursin("(", x)
        first_part = split(x, '(')[1]
        last_part = split(split(x, '(')[2], ')')[1]
        fla = Dict(
            string(atom.match) => parse.(Int, replace(split(first_part, r"[A-Za-z]+"), "" => "1")[2:end])[i] for (i, atom) in enumerate(eachmatch(r"[A-Za-z]+", first_part))
        )
        for (i, atom) in enumerate(eachmatch(r"[A-Za-z]+", last_part))
            if haskey(fla, atom.match)
                fla[atom.match] += 7 * parse.(Int, replace(split(last_part, r"[A-Za-z]+"), "" => "1")[2:end])[i]
            else
                fla[atom.match] = 7 * parse.(Int, replace(split(last_part, r"[A-Za-z]+"), "" => "1")[2:end])[i]
            end
        end
    else
        fla = Dict(
            string(atom.match) => parse.(Int, replace(split(x, r"[A-Za-z]+"), "" => "1")[2:end])[i] for (i, atom) in enumerate(eachmatch(r"[A-Za-z]+", x))
        )
    end
    return fla
end


"""
$(TYPEDSIGNATURES)

Add metabolic reactions to the model.
"""
function extend_model!(model, dfs)

    gs = String[]
    ms = RheaReactions.RheaMetabolite[]

    for df in dfs

        rid = parse(Int64, split(first(df.RHEA_ID), ':')[2])
        grr = String.(df.Protein[:])
        stoich = Int.(df.Stoichiometry[:])
        append!(gs, grr)

        if haskey(model.reactions, string(rid)) # isozyme

            push!(model.reactions[string(rid)].gene_association_dnf, grr)

        else # first time seeing this reaction
            rxn = get_reaction(rid)
            coeff_mets = get_reaction_metabolites(rid)
            stoichiometry = Dict(
                string(v.accession) => s
                for (s, v) in coeff_mets
            )

            append!(ms, last.(coeff_mets))

            ecs = isnothing(rxn.ec) ? df.EC : [rsplit(x, '/'; limit=2)[2] for x in rxn.ec]
            name = rxn.name

            #direction 
            reversibility_index_threshold = 5 
            rev_ind = ismissing(first(df.RevIndex)) ? nothing : first(df.RevIndex) 

            if isnothing(rev_ind) || (abs(rev_ind) <= reversibility_index_threshold)
                lb = -1000
                ub = 1000
            elseif rev_ind < -reversibility_index_threshold # forward
                lb = 0
                ub = 1000
            elseif rev_ind > reversibility_index_threshold # reverse
                lb = -1000
                ub = 0
            end


            model.reactions[string(rid)] = CM.Reaction(;
                name=name,
                lower_bound = lb,
                upper_bound = ub,
                stoichiometry = stoichiometry,
                gene_association_dnf = [grr],
                annotations = Dict(
                    "REACTION" => [rxn.equation],
                    "EC" => ecs,
                    "KEGG" => unique(df.KEGG_ID)
                ),
            )

        end
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

    # add genes
    for g in gs
        haskey(model.genes, g) && continue
        model.genes[g] = CM.Gene(; name=g)
    end

    model

end

function gapfill!(model)

    df = DataFrame(CSV.File("data/model/gapfilling_reactions.csv"))
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
                "reason" => ["gapfilling"]
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

"""
$(TYPEDSIGNATURES)

Get the SAPIG ID of a gene.
"""
function id_tag(gene)
    id_tag = Dict{String,String}()
    open("data/databases/ST398.txt") do io
        locus_tag = ""
        id = ""
        for ln in eachline(io)
            if startswith(ln, '>')
                id = split(ln; limit=2)[1][2:end]
                locus_tag = split(split(ln, "locus_tag=")[2], ']'; limit=2)[1]
            end
            id_tag[id] = locus_tag
        end
    end
    return id_tag[gene]
end
export id_tag

"""
$(TYPEDSIGNATURES)

Make model with ec number and gene ids as reaction names
"""
function change_reaction_names(model)
    id_tag = Dict{String,String}()
    open("data/databases/ST398.txt") do io
        locus_tag = ""
        id = ""
        for ln in eachline(io)
            if startswith(ln, '>')
                id = split(ln; limit=2)[1][2:end]
                locus_tag = split(split(ln, "locus_tag=")[2], ']'; limit=2)[1]
            end
            id_tag[id] = locus_tag
        end
    end
    tag_id = Dict(y => x for (x, y) in id_tag)

    eggnog = DataFrame(CSV.File("data/databases/MM_hla1w34o.emapper.annotations.tsv"))
    eggnog_dict = Dict(Pair.(eggnog.query, eggnog.Preferred_name))

    rhea_id_gene_id = Dict{String,String}()
    for (r, rxn) in model.reactions
        startswith(r,"EX") && continue
        g_name = ""
        if !isnothing(rxn.annotations) && haskey(rxn.annotations, "EC")
            for ec in split(rxn.annotations["EC"][1])
                g_name *= "$ec, "
            end
        end
        if !isnothing(rxn.gene_association_dnf)
            rxn.gene_association_dnf == [["g1"]] && continue
            for g in unique([eggnog_dict[tag_id[g]] for g in vcat(rxn.gene_association_dnf...)])
                g_name *= "$g, "
            end
        end
        rhea_id_gene_id[r] = g_name
    end

    g_model = deepcopy(model)
    for (rid, g_name) in rhea_id_gene_id
        length(g_name)<1 && continue
        g_model.reactions[rid].name = g_name
    end
    return g_model
end
export change_reaction_names

"""
$(TYPEDSIGNATURES)

Add source reactions to the model.
"""
function add_sources!(model)
    df = DataFrame(CSV.File("data/model/exchanges/sources.csv"))
    i = 0 
    for row in eachrow(df)
        i += 1
        mid = row.CHEBI
        chebi = split(mid,":")[2]
        if i ∈ [1,2,3,4,5,6,7,8,9]  
            mid *= "_e"
        end
        println(mid)
        model.reactions["EX_$chebi"] = CM.Reaction(
            ;
            name = "$(row.Name) exchange",
            lower_bound = 0.0,
            upper_bound = 1000.0,
            stoichiometry = Dict(mid => 1),
        )
        model.metabolites[mid] = deepcopy(model.metabolites[row.CHEBI])
        model.metabolites[mid].compartment = "external"
    end
    model
end

"""
$(TYPEDSIGNATURES)

Add source reactions to the model.
"""
function add_sinks!(model)
    df = DataFrame(CSV.File("data/model/exchanges/sinks.csv"))

    for row in eachrow(df)
        chebi = split(row.CHEBI,':')[2]
        model.reactions["EX_$chebi"] = CM.Reaction(
            ;
            name = "$(row.Name) exchange",
            lower_bound = -1000.0,
            upper_bound = 0.0,
            stoichiometry = Dict(row.CHEBI => 1),
        )
        model.metabolites[row.CHEBI] = deepcopy(model.metabolites[row.CHEBI])
        model.metabolites[row.CHEBI].compartment = "external"
    end
    model
end

function add_oxphos!(model)

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

function change_bounds!(model)
    df = DataFrame(CSV.File("data/model/unidirectional_reactions.csv"))

    for row in eachrow(df)
        model.reactions[string(row.RHEA_ID)].lower_bound = row.LOWER_BOUND 
        model.reactions[string(row.RHEA_ID)].upper_bound = row.UPPER_BOUND 
    end
    model
end
