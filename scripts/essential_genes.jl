using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using CairoMakie
using HiGHS, JSON, CSV
using JSONFBCModels, DataFrames, XLSX
using Latexify


model, reaction_isozymes = build_model()

model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_15903"].upper_bound = 10 #limit glucose

# ess_rxn = String[] 
# for (r,rxn) in model.reactions 
#     r == "biomass" && continue
#     lb = rxn.lower_bound 
#     ub = rxn.upper_bound 
#     model.reactions[r].lower_bound = 0 
#     model.reactions[r].upper_bound = 0 
#     sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
#     if isnothing(sol) || sol.objective < 1e-5 
#         push!(ess_rxn,r)
#     end
#     model.reactions[r].lower_bound = lb 
#     model.reactions[r].upper_bound = ub 
# end


# df = DataFrame(ID=String[],Name=String[],Stoichiometry=String[],Pathway=String[])
# for r in sort(ess_rxn)
#     println(r)
#     stoich = !haskey(model.reactions[r].annotations,"REACTION") ? "" : model.reactions[r].annotations["REACTION"][1]
#     fwrd = model.reactions[r].upper_bound > 0 && model.reactions[r].lower_bound >= 0
#     rvrs = model.reactions[r].upper_bound <= 0 && model.reactions[r].lower_bound < 0
#     bi = model.reactions[r].upper_bound > 0 && model.reactions[r].lower_bound < 0
#     if stoich == "" 
#         subs = [m for (m,s) in A.reaction_stoichiometry(model,r) if s<0]
#         prods = [m for (m,s) in A.reaction_stoichiometry(model,r) if s>0]
#         if fwrd || bi 
#             for m in subs 
#                 stoich *= A.metabolite_name(model,m)*" + "
#             end
#             stoich *= fwrd ? " => " : " <=> " 
#             for m in prods 
#                 stoich *= A.metabolite_name(model,m)*" + "
#             end
#         else 
#             for m in prods 
#                 stoich *= A.metabolite_name(model,m)*" + "
#             end
#             stoich *= fwrd ? " => " : " <=> " 
#             for m in subs 
#                 stoich *= A.metabolite_name(model,m)*" + "
#             end
#         end
#     elseif stoich != "" && rvrs 
#         stoich = string(split(stoich, " <=> "))[2]*" => "* string(split(stoich, " <=> "))[1]
#     end
#     push!(
#         df,
#         [
#             r,
#             isnothing(model.reactions[r].name) ? "" : model.reactions[r].name,
#             stoich,
#             !haskey(model.reactions[r].annotations,"Pathway") ? "" : join(model.reactions[r].annotations["Pathway"])
#         ]
#     )
# end

# df


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
                    "ID=$g" ∈ exp_df.Gene_ID ? "\\textcolor{red}{$g}" : g,
                    g == "" ? "" : g_name[g],
                    isnothing(model.reactions[r].name) ? r : model.reactions[r].name,
                    haskey(model.reactions[r].annotations,"Pathway") ? join(model.reactions[r].annotations["Pathway"]) : ""
                ]
            )
        end
    end
end

unique!(df)

new_gid = String[] 
for g in unique(df.GeneID)
    n_g = length([gid for gid in df.GeneID if gid==g])
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

[g for g in df.GeneID if g ∈ exp_df.Gene_ID]

newdf = filter(row -> row.GeneID ∈ exp_df.Gene_ID,df)

newdf = filter(row -> contains(row.Pathway,"mino acid"),df)
unique(newdf.GeneID)
