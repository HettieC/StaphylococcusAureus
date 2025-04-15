
function add_periplasm_transporters!(model)

    # all extracellular move with diffusion into periplasm
    all_exchange_metabolites =
        String.(DataFrame(
            CSV.File(
                joinpath(pkgdir(@__MODULE__), "data", "model", "exchanges", "sinks.csv"),
            ),
        ).CHEBI)
    append!(all_exchange_metabolites, String.(DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "model", "exchanges", "sources.csv"),
        ),
    ).CHEBI)
    )

    # bidirectional by default
    for mid in all_exchange_metabolites
        nm = A.metabolite_name(model, mid)

        model.metabolites[mid*"_p"] = deepcopy(model.metabolites[mid])

        model.metabolites[mid*"_p"].compartment = "Periplasm"

        model.reactions["DF_"*mid] = CM.Reaction(;
            name="Diffusion $nm",
            stoichiometry=Dict(mid * "_e" => -1, mid * "_p" => 1),
            lower_bound=-1000.0,
            upper_bound=1000,
            annotations=Dict("SBO" => ["SBO_0000284"]),
        )
        if !haskey(model.metabolites,mid * "_p")
            model.metabolites[mid * "_p"] = deepcopy(model.metabolites[mid])
            model.metabolites[mid * "_p"].compartment = "periplasm"
        end
        if !haskey(model.metabolites,mid * "_e")
            model.metabolites[mid * "_e"] = deepcopy(model.metabolites[mid])
            model.metabolites[mid * "_e"].compartment = "external"
        end
    end

end

function add_membrane_transporters!(model)

    df = DataFrame(
        CSV.File(joinpath(pkgdir(@__MODULE__), "data", "model", "transporters.csv")),
    )

    gs = String[]
    ms = String[]

    # abc transporters
    abcs = @subset(df, :Type .== "ABC")
    for g in groupby(abcs, [:CHEBI, :Subunit])
        all(x -> ismissing(x), g.Protein) && continue
        mid = first(g.CHEBI)
        mid != "CHEBI:15379" && continue

        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(filter(x -> !ismissing(x),g.Protein))
            append!(gs, iso)
            if all(x -> ismissing(x), g.Subunit)
                for gid in iso
                    add_abc!(model, mid, [gid], 1)
                end
            else
                ss = parse.(Float64, string.(g.Stoichiometry))
                add_abc!(model, mid, iso, ss)
            end
        else
            @warn "$mid not in model (ABC)"
        end
    end

    # PTS transporters
    pts = @subset(df, :Type .== "PTS")
    for g in groupby(pts, [:CHEBI, :Subunit])
        all(x -> ismissing(x), g.Protein) && continue
        mid = first(g.CHEBI)
        mid != "CHEBI:15379" && continue

        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(filter(x -> !ismissing(x),g.Protein))
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_pts!(model, mid, iso, ss)
        else
            @warn "$mid not in model (PTS)"
        end
    end

    # symport
    symport = @subset(df, :Type .== "Symport")
    for g in groupby(symport, [:CHEBI, :Subunit])
        all(x -> ismissing(x), g.Protein) && continue
        mid1, mid2 = sort(split(first(g.CHEBI), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(filter(x -> !ismissing(x),g.Protein))
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_symport!(model, mid1, mid2, iso, ss)
        else
            @warn "$mid1 or $mid2 not in model (symport)"
        end
    end

    # antiport
    antiport = @subset(df, :Type .== "Antiport")
    for g in groupby(antiport, [:CHEBI, :Subunit])
        all(x -> ismissing(x), g.Protein) && continue
        mid1, mid2 = sort(split(first(g.CHEBI), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(filter(x -> !ismissing(x),g.Protein))
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_antiport!(model, mid1, mid2, iso, ss)
        else
            @warn "$mid1 or $mid2 not in model (antiport)"
        end
    end

    # permease (the default as well)
    permease = @subset(df, :Type .== "Permease")
    for g in groupby(permease, [:CHEBI, :Subunit])
        all(x -> ismissing(x), g.Protein) && continue
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(filter(x -> !ismissing(x),g.Protein))
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_permease!(model, mid, iso, ss)
        else
            @warn "$mid not in model (permease)"
        end
    end

    # add default permeases - only reactions that did not get another transporter
    all_exchange_metabolites =
        String.(DataFrame(
            CSV.File(
                joinpath(pkgdir(@__MODULE__), "data", "model", "exchanges", "sinks.csv"),
            ),
        ).CHEBI)
    append!(all_exchange_metabolites, String.(DataFrame(
        CSV.File(
            joinpath(pkgdir(@__MODULE__), "data", "model", "exchanges", "sources.csv"),
        ),
    ).CHEBI)
    )
    # note: Pi, Na, and H will not get a permease here, due to them being involved in the other porters
    missing_transporters = setdiff(all_exchange_metabolites, unique(ms))
    for mid in missing_transporters
        if mid in A.metabolites(model)
            add_permease!(model, mid, nothing, [1.0])
        end
    end

    for g in gs
        haskey(model.genes, g) && continue
        model.genes[g] = CM.Gene(; name=g)
    end
    # no need to add metabolites, because they should all already be in the model
    @assert all(in.(ms, Ref(A.metabolites(model))))
end

function add_abc!(model, mid, iso, ss)
    rid = "ABC_$mid"
    # isoz =
    #     isnothing(iso) ? nothing :
    #     X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(iso) || push!(model.reactions[rid].gene_association_dnf, iso)
    else
        st = Dict(
            "CHEBI:30616" => -1, # atp
            "CHEBI:15377" => -1, # water
            "CHEBI:43474" => 1, # pi
            "CHEBI:456216" => 1, # adp
            "CHEBI:15378" => 1,  # h+ 
            mid * "_p" => -1.0,
        )
        st[mid] = get(st, mid, 0) + 1.0 # handle the case when phosphate is transported
        model.reactions[rid] = CM.Reaction(;
            name="Transport $(A.metabolite_name(model, String(mid))) ABC",
            stoichiometry=st,
            objective_coefficient=0.0,
            lower_bound=0,
            upper_bound=1000,
            gene_association_dnf=isnothing(iso) ? nothing : [iso],
            annotations=Dict("SBO" => ["SBO_0000284"]),
        )
        if !haskey(model.metabolites,mid * "_p")
            model.metabolites[mid * "_p"] = deepcopy(model.metabolites[mid])
            model.metabolites[mid * "_p"].compartment = "periplasm"
        end
    end
end

function add_pts!(model, mid, iso, ss)
    lu_phospho = Dict(
        "CHEBI:506227" => "CHEBI:57513", # n-acetyl-glucosamine -> N-acetyl-D-glucosamine 6-phosphate
        "CHEBI:15903" => "CHEBI:58247", # glucose -> glucose 6 phosphate
        "CHEBI:17992" => "CHEBI:57723", # sucrose -> sucrose 6 phosphate
        "CHEBI:16899" => "CHEBI:61381", # mannitol -> D-mannitol 1-phosphate
        "CHEBI:28645" => "CHEBI:57634", # β-D-fructose -> beta-D-fructose 6-phosphate
    )

    rid = "PTS_$mid"
    # isoz =
    #     isnothing(iso) ? nothing :
    #     X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(isoz) || push!(model.reactions[rid].gene_association_dnf, iso)
    else
        model.reactions[rid] = CM.Reaction(;
            name="Transport $(A.metabolite_name(model, String(mid))) PTS",
            stoichiometry=Dict(
                "CHEBI:58702" => -1.0, # pep
                "CHEBI:15361" => 1.0, # pyr
                mid * "_p" => -1.0,
                lu_phospho[mid] => 1.0, # cytosol phospho metabolite
            ),
            objective_coefficient=0.0,
            lower_bound=0,
            upper_bound=1000,
            gene_association_dnf=isnothing(iso) ? nothing : [iso],
            annotations=Dict("SBO" => ["SBO_0000284"]),
        )
        if !haskey(model.metabolites,mid * "_p")
            model.metabolites[mid * "_p"] = deepcopy(model.metabolites[mid])
            model.metabolites[mid * "_p"].compartment = "periplasm"
        end
    end
end

function add_symport!(model, mid1, mid2, iso, ss)
    rid = "SYM_$(mid1)_$mid2"
    # isoz =
    #     isnothing(iso) ? nothing :
    #     X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(isoz) || push!(model.reactions[rid].gene_association_dnf, iso)
    else
        model.reactions[rid] = CM.Reaction(;
            name="Symport $(A.metabolite_name(model, String(mid1)))::$(A.metabolite_name(model, String(mid2)))",
            stoichiometry=Dict(
                mid1 * "_p" => -1.0,
                mid1 => 1.0,
                mid2 * "_p" => -1.0,
                mid2 => 1.0,
            ),
            objective_coefficient=0.0,
            lower_bound=0,
            upper_bound=1000,
            gene_association_dnf=isnothing(iso) ? nothing : [iso],
            annotations=Dict("SBO" => ["SBO_0000284"]),
        )
        if !haskey(model.metabolites,mid1 * "_p")
            model.metabolites[mid1 * "_p"] = deepcopy(model.metabolites[mid1])
            model.metabolites[mid1 * "_p"].compartment = "periplasm"
        end
        if !haskey(model.metabolites,mid2 * "_p")
            model.metabolites[mid2 * "_p"] = deepcopy(model.metabolites[mid2])
            model.metabolites[mid2 * "_p"].compartment = "periplasm"
        end
    end
end

function add_antiport!(model, mid1, mid2, iso, ss)
    rid = "ANTI_$(mid1)_$mid2"
    # isoz =
    #     isnothing(iso) ? nothing :
    #     X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(iso) || push!(model.reactions[rid].gene_association_dnf, iso)
    else
        model.reactions[rid] = CM.Reaction(;
            name="Antiport $(A.metabolite_name(model, String(mid1)))::$(A.metabolite_name(model, String(mid2)))",
            stoichiometry=Dict(
                mid1 * "_p" => 1.0,
                mid1 => -1.0,
                mid2 * "_p" => -1.0,
                mid2 => 1.0,
            ),
            objective_coefficient=0.0,
            lower_bound=0,
            upper_bound=1000,
            gene_association_dnf=isnothing(iso) ? nothing : [iso],
            annotations=Dict("SBO" => ["SBO_0000284"]),
        )
        if !haskey(model.metabolites,mid1 * "_p")
            model.metabolites[mid1 * "_p"] = deepcopy(model.metabolites[mid1])
            model.metabolites[mid1 * "_p"].compartment = "periplasm"
        end
        if !haskey(model.metabolites,mid2 * "_p")
            model.metabolites[mid2 * "_p"] = deepcopy(model.metabolites[mid2])
            model.metabolites[mid2 * "_p"].compartment = "periplasm"
        end
    end
end

function add_permease!(model, mid, iso, ss)
    rid = "PERM_$mid"
    # isoz =
    #     isnothing(iso) ? nothing :
    #     X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        isnothing(isoz) || push!(model.reactions[rid].gene_association_dnf, iso)
    else
        model.reactions[rid] = CM.Reaction(;
            name="Permease $(A.metabolite_name(model, String(mid)))",
            stoichiometry=Dict(mid * "_p" => -1.0, mid => 1.0),
            objective_coefficient=0.0,
            lower_bound=-1000,
            upper_bound=1000,
            gene_association_dnf=isnothing(iso) ? nothing : [iso],
            annotations=Dict("SBO" => ["SBO_0000284"]),
        )
        if !haskey(model.metabolites,mid* "_p")
            model.metabolites[mid * "_p"] = deepcopy(model.metabolites[mid])
            model.metabolites[mid * "_p"].compartment = "periplasm"
        end
    end
end
