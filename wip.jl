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
tag_id = Dict(y=>x for (x,y) in id_tag)

eggnog = DataFrame(CSV.File("data/databases/MM_hla1w34o.emapper.annotations.tsv"))
eggnog_dict = Dict(Pair.(eggnog.query,eggnog.Preferred_name))

rhea_id_gene_id = Dict{String,String}()
for (r,rxn) in model.reactions 
    isnothing(rxn.gene_association_dnf) && continue 
    g_name = ""
    if !isnothing(rxn.annotations) && haskey(rxn.annotations,"EC")
        for ec in split(rxn.annotations["EC"][1])
            g_name *= "$ec, "
        end
    end
    for g in unique([eggnog_dict[tag_id[g]] for g in vcat(rxn.gene_association_dnf...)])
        g_name *= "$g, "
    end
    rhea_id_gene_id[r] = g_name
end

g_model=deepcopy(model)
for (rid,g_name) in rhea_id_gene_id
    g_model.reactions[rid].name = g_name 
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


get_reactions_with_ec("2.7.2.8")

get_reaction(14629)


id_tag["lcl|AM990992.1_prot_CAQ48634.1_194"]




