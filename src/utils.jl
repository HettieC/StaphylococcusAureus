"""
$(TYPEDSIGNATURES)

Parse the chemical formula from a RheaReaction.formula
"""
function parse_formula(x::Union{Nothing,String})
    isnothing(x) && return nothing
    x == "" && return nothing

    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"
    for m in eachmatch(pattern, x)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end
    return res
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
            # only use CHEBI metabolites
            any(((k,v),)->!startswith(k,"CHEBI"),stoichiometry) && continue

            append!(ms, last.(coeff_mets))

            if isnothing(rxn)
                println(rid)
                continue 
            end
            ecs = isnothing(rxn.ec) ? df.EC : [rsplit(x, '/'; limit=2)[2] for x in rxn.ec]
            name = rxn.name

            model.reactions[string(rid)] = CM.Reaction(;
                name=name,
                lower_bound = -1000,
                upper_bound = 1000,
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
            "Cytosol",
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

    return model

end

function gapfill!(model)

    df = DataFrame(CSV.File("data/model/reactions/gapfilling_reactions.csv"))
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
            gene_association_dnf = [["g1"]]
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
        if !isnothing(rxn.gene_association_dnf)
            rxn.gene_association_dnf == [["g1"]] && continue
            for g in unique([eggnog_dict[tag_id[g]] for g in vcat(rxn.gene_association_dnf...)])
                g == "-" && continue
                if g_name == "" 
                    g_name *= "$g"
                else 
                    g_name *=", $g"
                end
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
    for row in eachrow(df)
        mid = row.CHEBI
        chebi = split(mid,":")[2]
        mid *= "_e"
        model.reactions["EX_$chebi"] = CM.Reaction(
            ;
            name = "$(row.Name) exchange",
            lower_bound = 0,
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
        mid = row.CHEBI
        chebi = split(row.CHEBI,':')[2]
        mid *= "_e"
        if haskey(model.reactions,"EX_$chebi")
            model.reactions["EX_$chebi"].lower_bound = -1000 
        else
            model.reactions["EX_$chebi"] = CM.Reaction(
                ;
                name = "$(row.Name) exchange",
                lower_bound = -1000.0,
                upper_bound = 0.0,
                stoichiometry = Dict(mid => 1),
            )
            model.metabolites[mid] = deepcopy(model.metabolites[row.CHEBI])
            model.metabolites[mid].compartment = "external"
        end
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
    df = DataFrame(CSV.File("data/model/reactions/unidirectional_reactions.csv"))

    for row in eachrow(df)
        !haskey(model.reactions,string(row.RHEA_ID)) && continue
        model.reactions[string(row.RHEA_ID)].lower_bound = row.LOWER_BOUND 
        model.reactions[string(row.RHEA_ID)].upper_bound = row.UPPER_BOUND 
    end
    model
end

function get_gene_product_molar_mass(gids)
    AA_mass = Dict(
        'A' => 71.0788,
        'R' => 156.1875,
        'N' => 114.1038,
        'D' => 115.0886,
        'C' => 103.1388,
        'Q' => 128.1307,
        'E' => 129.1155,
        'G' => 57.0519,
        'H' => 137.1411,
        'I' => 113.1594,
        'L' => 113.1594,
        'K' => 128.1741,
        'M' => 131.1926,
        'F' => 147.1766,
        'P' => 97.1167,
        'S' => 87.0782,
        'T' => 101.1051,
        'W' => 186.2132,
        'Y' => 163.1760,
        'V' => 99.1326,
    )
    seqs = Dict{String,String}()
    open("data/sequence.txt","r") do io
        prev = ""
        for ln in eachline(io)
            if startswith(ln,'>')
                seqs[ln[2:end]] = ""
                prev = ln[2:end]
            else 
                seqs[prev] = seqs[prev]*ln 
            end
        end
    end

    gene_product_molar_mass = Dict{String,Float64}()
    for gid in gids 
        gene_product_molar_mass[gid] = sum([AA_mass[aa] for aa in seqs[gid]]) / 1000 # convert to kDa
    end

    return gene_product_molar_mass
end

export get_gene_product_molar_mass

function add_isozymes!(reaction_isozymes,kcat_dict,dfs)
    for df in dfs
        rid = split(first(df.RHEA_ID),':')[2]
        !haskey(kcat_dict, "$(rid)_f") && continue # skip if no kcat data available

        grr = String.(df.Protein[:])
        stoich = Int.(df.Stoichiometry[:])
        
        d = get!(reaction_isozymes, rid, Dict{String,Isozyme}())
        i = length(reaction_isozymes[rid])
        d["isozyme_"*string(i+1)] = Isozyme( # each isozyme gets a unique name
            gene_product_stoichiometry=Dict(grr .=> stoich),
            kcat_forward=maximum([kcat_dict["$(rid)_f"][g] for g in grr]) * 3.6, #mmol/h
            kcat_reverse=maximum([kcat_dict["$(rid)_r"][g] for g in grr]) * 3.6,
        )
    end
    return reaction_isozymes
end

function get_reaction_isozymes()
    turnup_df = DataFrame(XLSX.readtable("data/turnup/output1.xlsx", "Sheet1"))
    turnup_df = vcat(turnup_df, DataFrame(XLSX.readtable("data/turnup/output2.xlsx", "Sheet1")))
    turnup_df = vcat(turnup_df, DataFrame(XLSX.readtable("data/turnup/output3.xlsx", "Sheet1")))
    turnup_df = vcat(turnup_df, DataFrame(XLSX.readtable("data/turnup/output4.xlsx", "Sheet1")))
    turnup_df = vcat(turnup_df, DataFrame(XLSX.readtable("data/turnup/output5.xlsx", "Sheet1")))
    turnup_df = vcat(turnup_df, DataFrame(XLSX.readtable("data/turnup/output6.xlsx", "Sheet1")))
    turnup_df = vcat(turnup_df, DataFrame(XLSX.readtable("data/turnup/output7.xlsx", "Sheet1")))

    rxn_seq_subs_prods = DataFrame(CSV.File("data/model/isozymes/reaction_sequence_subs_prods.csv"))
    
    turnup_df = DataFrames.rename!(turnup_df,
        "kcat [s^(-1)]" => "kcat",
    )

    kcat_df = insertcols(rxn_seq_subs_prods, 5, :kcat => turnup_df.kcat)
    kcat_dict = Dict{String,Dict{String,Float64}}()
    for row in eachrow(kcat_df)
        if !ismissing(row.kcat)
            if !haskey(kcat_dict, row.rxn)
                # put into 1/h instead of 1/s by multiplying by 3600
                kcat_dict[row.rxn] = Dict(row.locus_tag => row.kcat * 3.6)
            else
                kcat_dict[row.rxn][row.locus_tag] = row.kcat * 3.6
            end
        end
    end
    
    df = DataFrame(CSV.File("data/model/reactions/metabolic_reactions.csv"))
    
    heteros = @rsubset(df, !iszero(:Isozyme))
    @select!(heteros, :RHEA_ID, :Protein, :Stoichiometry, :Isozyme)
    
    homos = @rsubset(df, iszero(:Isozyme))
    @select!(homos, :RHEA_ID, :Protein, :Stoichiometry)
    
    ghomos = groupby(homos, [:RHEA_ID, :Protein])
    gheteros = groupby(heteros, [:RHEA_ID, :Isozyme])
    
    reaction_isozymes = Dict{String,Dict{String,Isozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
    add_isozymes!(reaction_isozymes,kcat_dict,ghomos)
    add_isozymes!(reaction_isozymes,kcat_dict,gheteros)
    
    
    return reaction_isozymes, kcat_dict
end
export get_reaction_isozymes

function add_fake_isozymes!(model,reaction_isozymes)
    # use a fake gene g1 for all metabolic reactions with not grr

    avg_kcat = sum(vcat([b.kcat_forward for (x,y) in reaction_isozymes for (a,b) in y]...))/length([b.kcat_forward for (x,y) in reaction_isozymes for (a,b) in y])

    for (r,rxn) in model.reactions 
        haskey(reaction_isozymes,r) && continue 
        startswith(r,"EX") && continue 
        startswith(r,"DF") && continue 
        r == "biomass" && continue
        if rxn.gene_association_dnf == [["g1"]]
            reaction_isozymes[r] = Dict("isozyme_1" => Isozyme(
                gene_product_stoichiometry = Dict("g1" => 1),
                kcat_forward = startswith(r,"PERM") ? 60 : avg_kcat,
                kcat_reverse = startswith(r,"PERM") ? 60 : avg_kcat
            )
            )
        else
            println("fix isozymes of $r")
        end
    end

    oxphos_reactions = ["Ndh2", "Sdh", "Mqo", "Lqo", "cyt_aa3", "cyt_bd", "cyt_bo3"]
    for rid in oxphos_reactions
        reaction_isozymes[rid] = Dict("isoyzme_1" => Isozyme(
            gene_product_stoichiometry = Dict("g1" => 1), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = avg_kcat,
            kcat_reverse = avg_kcat,
            )
        )
    end
    return reaction_isozymes
end

function find_ec_info(ec_number::String)
    df = DataFrame(CSV.File("data/databases/uniprot/s_aureus.csv"))
    a = @subset df @byrow begin
        !ismissing(:e_c_number) && contains(:e_c_number, ec_number)
    end
    return a
end
export find_ec_info

function add_special_isozymes!(reaction_isozymes,kcat_dict,model)
    for (rid,rxn) in model.reactions 
        if !haskey(reaction_isozymes,rid) && haskey(kcat_dict,rid*"_f")
            reaction_isozymes[rid] = Dict(
                "isozyme_"*string(i) => Isozyme(
                    gene_product_stoichiometry = Dict(g => 1 for g in grr),
                    kcat_forward = maximum([kcat_dict["$(rid)_f"][g] for g in grr]) * 3.6,
                    kcat_reverse = maximum([kcat_dict["$(rid)_r"][g] for g in grr]) * 3.6
                )
                for (i,grr) in enumerate(rxn.gene_association_dnf)
            )
        end
    end
    #make ATPS faster
    reaction_isozymes["ATPS"] = Dict(
        id => Isozyme(
            iso.gene_product_stoichiometry,
            iso.kcat_forward * 10,
            iso.kcat_reverse
        )
        for (id,iso) in reaction_isozymes["ATPS"]
    )
    return reaction_isozymes
end

function add_genes!(model)
    for (r,rxn) in model.reactions
        isnothing(rxn.gene_association_dnf) && continue 
        for gid in vcat(rxn.gene_association_dnf...)
            if !haskey(model.genes,gid)
                model.genes[gid] = CM.Gene(;name = gid)
            end
        end
    end
    return model 
end

function _add_kegg_info!(model)
    for (r,rxn) in model.reactions 
        println(r)
        if !isnothing(rxn.name) && !isnothing(tryparse(Int,r))
            rxn.annotations["NameDB"] = ["Rhea"]
            kegg_info = get_kegg_info(rxn.annotations["KEGG"][1])
            model.reactions[r].annotations["Pathway"] = isnothing(kegg_info["pathway"]) ? [""] : kegg_info["pathway"]
        elseif isnothing(tryparse(Int,r))
            rxn.annotations["NameDB"] = ["none"]
        elseif haskey(rxn.annotations,"KEGG") && rxn.annotations["KEGG"] != ["R"]
            println(r)
            rxn.annotations["NameDB"] = ["KEGG"]
            kegg_info = get_kegg_info(rxn.annotations["KEGG"][1])
            isnothing(kegg_info) && continue
            model.reactions[r].name = kegg_info["name"]
            model.reactions[r].annotations["Pathway"] = isnothing(kegg_info["pathway"]) ? [""] : kegg_info["pathway"]
        end
    end
    return model 
end


