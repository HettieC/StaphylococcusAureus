using CSV, DataFrames, DataFramesMeta
using COBREXA, AbstractFBCModels
import AbstractFBCModels.CanonicalModel as CM
import COBREXA as X
using DocStringExtensions
using RheaReactions
using JSONFBCModels

include("src/utils.jl")
include("src/reconstruct.jl")

model = build_model()

save_model(convert(JSONFBCModels.JSONFBCModel,model),"data/model.json")

# make model with gene ids as reaction names
df = DataFrame(CSV.File("data/model/gene_ids.csv"))
g_model = deepcopy(model)
for row in eachrow(df)
    g_model.reactions[string(row.RHEA_ID)].name = row.Gene_ID 
end

save_model(convert(JSONFBCModels.JSONFBCModel,g_model),"escher_model.json")


####################################
id_tag = Dict{String,String}()
open("data/databases/ST398.txt") do io
    locus_tag = ""
    id = ""
    for ln in eachline(io)
        if startswith(ln, '>')
            id = split(ln;limit=2)[1][2:end]
            locus_tag = split(split(ln, "locus_tag=")[2], ']'; limit = 2)[1]
        end
        id_tag[id] = locus_tag
    end
end



function add_membrane_transporters!(model)
   
    df = DataFrame(
        CSV.File(
            joinpath("data", "model", "transporters.csv"),
        ),
    )

    gs = String[]
    ms = String[]

    # abc transporters
    abcs = @subset(df, :Type .== "ABC")
    for g in groupby(abcs, [:CHEBI, :Subunit])
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
            ss = parse.(Float64, string.(g.Stoichiometry))
            append!(gs, iso)
            add_abc!(model, mid, iso, ss)
        else
            @warn "$mid not in model (ABC)"
        end
    end

    # PTS transporters
    pts = @subset(df, :Type .== "PTS")
    for g in groupby(pts, [:CHEBI, :Subunit])
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
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
        mid1, mid2 = sort(split(first(g.CHEBI), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(g.Protein)
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
        mid1, mid2 = sort(split(first(g.CHEBI), "/")) # to make rid unique
        if mid1 in A.metabolites(model) && mid2 in A.metabolites(model)
            push!(ms, mid1)
            push!(ms, mid2)
            iso = string.(g.Protein)
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
        mid = first(g.CHEBI)
        if mid in A.metabolites(model)
            push!(ms, mid)
            iso = string.(g.Protein)
            append!(gs, iso)
            ss = parse.(Float64, string.(g.Stoichiometry))
            add_permease!(model, mid, iso, ss)
        else
            @warn "$mid not in model (permease)"
        end
    end

    # add default permeases
    all_exchange_metabolites =  DataFrame(
        CSV.File(
            joinpath("data", "model", "exchange_metabolites.csv"),
        ),
    ).CHEBI
    # note: Pi, Na, and H will not get a permease here, due to them being involved in the other porters
    missing_transporters = setdiff(all_exchange_metabolites, unique(ms))
    for mid in missing_transporters
        if mid in A.metabolites(model)
            add_permease!(model, mid, ["Missing"], [1.0])
        end
    end

    add_genes!(model, gs)
    # no need to add metabolites, because they should all already be in the model
    @assert all(in.(ms, Ref(A.metabolites(model))))
end

function add_abc!(model, mid, iso, ss)
    rid = "ABC_$mid"
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $(A.metabolite_name(model, String(mid))) ABC",
            stoichiometry = Dict(
                "CHEBI:30616" => -1, # atp
                "CHEBI:15377" => -1, # water
                "CHEBI:43474" => 1, # pi
                "CHEBI:456216" => 1, # adp
                "CHEBI:15378" => 1,  # h+ 
                mid * "_p" => -1.0,
                mid => 1.0, # cytosol
            ),
            objective_coefficient = 0.0,
            lower_bound = 0,
            upper_bound = 1000,
            gene_association = [isoz],
        )
    end
end

function add_pts!(model, mid, iso, ss)
    lu_phospho = Dict(
        "CHEBI:506227" => "CHEBI:57513", # n-acetyl-glucosamine -> N-acetyl-D-glucosamine 6-phosphate
        "CHEBI:15903" => "CHEBI:58247", # glucose -> glucose 6 phosphate
        "CHEBI:17992" => "", # sucrose -> 
        "CHEBI:16899" => "CHEBI:61381", # mannitol -> D-mannitol 1-phosphate
        "CHEBI:17997" => "", # Nitrogen -> 

    )

    rid = "PTS_$mid"
    isoz = X.Isozyme(; gene_product_stoichiometry = Dict(iso .=> ss))
    if haskey(model.reactions, rid)
        push!(model.reactions[rid].gene_association, isoz)
    else
        model.reactions[rid] = Reaction(
            name = "Transport $(A.metabolite_name(model, String(mid))) PTS",
            stoichiometry = Dict(
                "CHEBI:58702" => -1.0, # pep
                "CHEBI:15361" => 1.0, # pyr
                mid * "_p" => -1.0,
                lu_phospho[mid] => 1.0, # cytosol phospho metabolite
            ),
            objective_coefficient = 0.0,
            lower_bound = 0,
            upper_bound = 1000,
            gene_association = [isoz],
        )
    end
end






