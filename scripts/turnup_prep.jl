using DataFrames, CSV, StaphylococcusAureus, XLSX

model = build_model()

chebi_df = DataFrame(CHEBI=String[],name = String[], smiles=String[],inchi=String[],inchikey=String[])
open("data/databases/chebi/chebi_core.obo","r") do io 
    i = 0
    chebi = "" 
    Name = ""
    smile = ""
    inchi_key = "" 
    inchi_val = ""
    for ln in eachline(io)
        i += 1
        i < 20 && continue
        if startswith(ln,"id: ")
            chebi = string(split(ln,"id: ")[2])
        elseif startswith(ln,"name: ")
            Name = string(split(ln; limit = 2)[2])
        elseif startswith(ln,"property_value: http://purl.obolibrary.org/obo/chebi/smiles ")
            smile = string(split(ln)[3][2:end-1])
        elseif startswith(ln, "property_value: http://purl.obolibrary.org/obo/chebi/inchi ")
            inchi_val = string(split(ln)[3][2:end-1])
        elseif startswith(ln, "property_value: http://purl.obolibrary.org/obo/chebi/inchikey ")
            inchi_key = string(split(ln)[3][2:end-1])
        elseif ln == "[Term]"
            push!(chebi_df, [chebi, Name, smile, inchi_val, inchi_key])
        end
    end
end
chebi_inchi_dict = Dict(Pair.(chebi_df.CHEBI,chebi_df.inchi))

seq_dict = Dict{String,String}()
id_tag = Dict{String,String}()
open("data/databases/ST398.txt") do io
    locus_tag = ""
    seq = ""
    id = ""
    for ln in eachline(io)
        if startswith(ln, '>')
            id = split(ln;limit=2)[1][2:end]
            locus_tag = split(split(ln, "locus_tag=")[2], ']'; limit = 2)[1]
            if !haskey(seq_dict, locus_tag)
                seq_dict[locus_tag] = ""
            end
        else
            seq = seq_dict[locus_tag]
            seq_dict[locus_tag] = "$seq$ln"
            id_tag[id] = locus_tag
        end
    end
end


df = DataFrame(
    rxn = String[],
    locus_tag = String[],
    Enzyme = String[],
    Substrates = String[],
    Products = String[],
)
for (id, r) in model.reactions
    if !isnothing(r.gene_association_dnf) && r.gene_association_dnf != [["g1"]]
        subs = ""
        prods = ""
        for (x, y) in r.stoichiometry
            if contains(x,"_")
                x = split(x,'_')[1]
            end
            !haskey(chebi_inchi_dict,x) && continue
            if y < 0
                subs = "$(chebi_inchi_dict[x]);$subs"
            else
                prods = "$(chebi_inchi_dict[x]);$prods"
            end
        end

        tags = vcat(r.gene_association_dnf...)

        for t in tags
            push!(df, ("$(id)_f", t, seq_dict[t], subs, prods))
            push!(df, ("$(id)_r", t, seq_dict[t], prods, subs))
        end
    end
end

CSV.write("data/model/isozymes/reaction_sequence_subs_prods.csv",df)
select!(df, [:Enzyme, :Substrates, :Products])

CSV.write("data/turnup/sequence_subs_prods.csv", df; delim = ',')
XLSX.writetable("data/turnup/turnup_input1.xlsx",df[1:499,:])
XLSX.writetable("data/turnup/turnup_input2.xlsx",df[500:999,:])
XLSX.writetable("data/turnup/turnup_input3.xlsx",df[1000:1499,:])
XLSX.writetable("data/turnup/turnup_input4.xlsx",df[1499:end,:])

unique(turnup_df)
unique(rxn_seq_subs_prods)
DataFrames.rename!(turnup_df,
        [
            "substrates" => "Substrates",
            "enzyme" => "Enzyme",
            "products" => "Products"
        ]
    )
innerjoin(rxn_seq_subs_prods,turnup_df, on = [:Enzyme,:Substrates,:Products])
