using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using HiGHS, JSON, CSV
using JSONFBCModels, DataFrames, XLSX
using Latexify

model, reaction_isozymes = build_model()

model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_15903"].upper_bound = 10 #limit glucose


# find essential genes on minimal medium
ko_dict = gene_knockouts(model, optimizer = HiGHS.Optimizer)

#### use richer medium 
# 
rich_medium = DataFrame(CSV.File("data/model/exchanges/rich_medium.csv"))
for row in eachrow(rich_medium)
    if haskey(model.reactions, "EX_$(split(row.CHEBI,":")[2])") 
        model.reactions["EX_$(split(row.CHEBI,":")[2])"].upper_bound = 1000
    else
        model.reactions["EX_rich_$(split(row.CHEBI,":")[2])"] = CM.Reaction(
            name = "Exchange of $(row.Name)",
            lower_bound = 0,
            upper_bound = 1000,
            stoichiometry = Dict(row.CHEBI => 1.0),
        )
    end
end
ko_dict = gene_knockouts(model, optimizer = HiGHS.Optimizer)

kos = [k for (k,v) in ko_dict if k != "g1" && (isnothing(v) || abs(v) < 1e-4)]


# get gene names
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
g_name = Dict(g => eggnog_dict[tag_id[g]] for g in A.genes(model) if haskey(tag_id,g))


# make df of all genes and reactions for supplementary
df = DataFrame(GeneID=String[],Name=String[],ReactionID=String[],Reaction=String[])
for (g,ko) in ko_dict
    g == "g1" && continue
    if abs(ko)<1e-4
        rxns = [r for (r,rxn) in model.reactions if !isnothing(rxn.gene_association_dnf) && g ∈ vcat(rxn.gene_association_dnf...)]
        for r in rxns
            push!(
                df,
                [
                    g,
                    g == "" ? "" : g_name[g],
                    r,
                    isnothing(model.reactions[r].name) ? r : model.reactions[r].name,
                ]
            )
        end
    end
end
df
unique!(df)
new_gid = String[] 
i = 0
for g in unique(df.GeneID)
    i += 1
    n_g = length([gid for gid in df.GeneID if gid==g])
    if n_g == 1 
        if iseven(i)
            push!(new_gid,"\\rowcolor[gray]{0.9}$g")
        else 
            push!(new_gid,g)
        end
    else 
        if iseven(i)
            push!(new_gid, "\\rowcolor[gray]{0.9}$g")
            append!(new_gid,repeat(["\\rowcolor[gray]{0.9}"],n_g-1))
        else
            push!(new_gid, g)
            append!(new_gid,repeat([""],n_g-1))
        end
    end
end
df.GeneID = new_gid
new_gname = String[]
for g in unique(df.Name)
    n_g = length([gid for gid in df.Name if gid==g])
    if n_g == 1 
        push!(new_gname,g)
    else 
        push!(new_gname, g)
        append!(new_gname,repeat([""],n_g-1))
    end
end
df.Name = new_gname
latexify(df; env = :table, booktabs = true, latex = false) |> print



#experimental data
exp_df = DataFrame(XLSX.readtable("data/experimental/pone.0089018.s006.xlsx", "Zero transposon inserts"))
rename!(exp_df, Symbol.(replace.(string.(names(exp_df)), Ref(r"\s"=>"_"))))
filter!(row -> string(split(row.Gene_ID,"=")[2]) ∈ A.genes(model),exp_df)
exp_df.Gene_ID = string.([split(g,"=")[2] for g in exp_df.Gene_ID])

nader_df = DataFrame(XLSX.readtable("data/experimental/Essential_gene_Nader.xlsx", "essential genes"))

unique(vcat(exp_df.Gene_ID, filter(g -> g ∈ A.genes(model), nader_df.locus_tag)))




# only experimental matches
exp_match = filter(row->row.GeneID ∈ exp_df.Gene_ID || row.GeneID ∈ nader_df.locus_tag, df)

filter!(row -> row.GeneID ∈ ["SAPIG1311","SAPIG1427","SAPIG2388"], exp_match)

new_gid = String[] 
for g in unique(exp_match.GeneID)
    n_g = length([gid for gid in exp_match.GeneID if gid==g])
    if n_g == 1 
        push!(new_gid,g)
    else 
        push!(new_gid, "\\multirow{$(n_g)}{=}{$g}")
        append!(new_gid,repeat([""],n_g-1))
    end
end
exp_match.GeneID = new_gid

new_name = String[] 
for g in unique(exp_match.Name)
    n_g = length([gid for gid in exp_match.Name if gid==g])
    if n_g == 1 
        push!(new_name,g)
    else 
        push!(new_name, "\\multirow{$(n_g)}{=}{$g}")
        append!(new_name,repeat([""],n_g-1))
    end
end
exp_match.Name = new_name

latexify(exp_match; env = :table, booktabs = true, latex = false) |> print
