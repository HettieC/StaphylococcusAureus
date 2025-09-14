"""
Add the KEGG reaction names and pathways. 
"""
function add_names_pathways!(model)
    dic = Dict(x => y for (x,y) in JSON.parsefile("data/model/reactions/kegg_names_pathways.json"))
    for (x,y) in dic
        !haskey(model.reactions,x) && continue
        model.reactions[x].name = y[1] 
        model.reactions[x].annotations["Pathway"] = y[2] == [""] ? [""] : [string(split(p,"  ")[2])*"; " for p in y[2]]
    end
end

function metanetx_annotate!(model)
    id_map = JSON.parsefile("data/databases/metanetx/id-mapper.json")
    for (k,v) in id_map 
        !haskey(v,"xrefs") && continue
        model.reactions[k].annotations["BiGG"] = unique([split(id,":")[2] for id in v["xrefs"] if startswith(id,"bigg")])
        model.reactions[k].annotations["metacyc"] = unique([split(id,":")[2] for id in v["xrefs"] if startswith(id,"metacyc")])
        model.reactions[k].annotations["seed"] = unique([split(id,":")[2] for id in v["xrefs"] if startswith(id,"seed")])
        model.reactions[k].annotations["sabiork"] = unique([split(id,":")[2] for id in v["xrefs"] if startswith(id,"sabiork")])
        model.reactions[k].annotations["RHEA"] = unique([split(id,":")[2] for id in v["xrefs"] if startswith(id,"rhea")])
        append!(model.reactions[k].annotations["KEGG"], unique([split(id,":")[2] for id in v["xrefs"] if startswith(id,"kegg")]))
    end
    return model 
end
