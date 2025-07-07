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


# find essential genes
ko_dict = gene_knockouts(model, optimizer = HiGHS.Optimizer)
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

#experimental data
exp_df = DataFrame(XLSX.readtable("data/experimental/pone.0089018.s006.xlsx", "Zero transposon inserts"))
rename!(exp_df, Symbol.(replace.(string.(names(exp_df)), Ref(r"\s"=>"_"))))
filter!(row -> string(split(row.Gene_ID,"=")[2]) ∈ A.genes(model),exp_df)
exp_df.Gene_ID = string.([split(g,"=")[2] for g in exp_df.Gene_ID])

nader_df = DataFrame(XLSX.readtable("data/experimental/Essential_gene_Nader.xlsx", "essential genes"))

# make df of gene and reaction
df = DataFrame(GeneID=String[],Name=String[],Reaction=String[],Pathway=String[])
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
                    isnothing(model.reactions[r].name) ? r : model.reactions[r].name,
                    haskey(model.reactions[r].annotations,"Pathway") ? join(model.reactions[r].annotations["Pathway"]) : ""
                ]
            )
        end
    end
end
df
unique!(df)


exp_match = filter(row->row.GeneID ∈ exp_df.Gene_ID,df)# || g ∈ nader_df.locus_tag])

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

new_pway = String[] 
for p in exp_match.Pathway 
    p = replace(p, "; Biosynthesis of secondary metabolites;" => "; ")
    p = replace(p, "Metabolic pathways; " => "")
    p = replace(p, "Fatty acid biosynthesis; Fatty acid metabolism" => "Fatty acid biosynthesis; ")
    p = replace(p, "Fatty acid degradation; Fatty acid metabolism" => "Fatty acid biosynthesis; ")
    p = replace(p, "; Biosynthesis of secondary metabolites;" => "; ")
    p = replace(p, "Histidine metabolism;  Biosynthesis of amino acids; " => "Histidine metabolism; ")
    push!(new_pway,p)
end
exp_match.Pathway = new_pway

latexify(exp_match; env = :table, booktabs = true, latex = false) |> print




new_gid = String[] 
for g in unique(df.GeneID)
    n_g = length([gid for gid in df.GeneID if gid==g])
    if g ∈ exp_df.Gene_ID || g ∈ nader_df.locus_tag
        g = "\\textcolor{red}{$g}"
    end
    if n_g == 1 
        push!(new_gid,g)
    else 
        push!(new_gid, "\\multirow{$(n_g)}{=}{$g}")
        append!(new_gid,repeat([""],n_g-1))
    end
end
df.GeneID = new_gid


new_pway = String[] 
for p in df.Pathway 
    p = replace(p, "; Biosynthesis of secondary metabolites;" => "; ")
    p = replace(p, "Metabolic pathways; " => "")
    p = replace(p, "Fatty acid biosynthesis; Fatty acid metabolism" => "Fatty acid biosynthesis; ")
    p = replace(p, "Fatty acid degradation; Fatty acid metabolism" => "Fatty acid biosynthesis; ")
    p = replace(p, "; Biosynthesis of secondary metabolites;" => "; ")
    p = replace(p, "Histidine metabolism;  Biosynthesis of amino acids; " => "Histidine metabolism; ")
    push!(new_pway,p)
end
df.Pathway = new_pway

latexify(df; env = :table, booktabs = true, latex = false) |> print

newdf = filter(row -> row.GeneID ∈ exp_df.Gene_ID,df)

newdf = filter(row -> contains(row.Pathway,"mino acid"),df)
unique(newdf.GeneID)
