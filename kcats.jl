using CSV, DataFrames, DataFramesMeta
using JSON

#### make isozyme df 
iso = DataFrame(CSV.File("data/isozymes.csv"))
rxns = DataFrame(CSV.File("data/model/metabolic_reactions.csv"))
select!(rxns, :KEGG_ID, :RHEA_ID, :Protein)

df = leftjoin(rxns, iso, on=[:Protein, :KEGG_ID])
CSV.write("reaction_isozymes.csv", df)

# ### get reaction turnover numbers kcat

# ### get sequence, substrates, products file
# seq_dict = Dict{String,String}()
# open("data/databases/ST398.txt") do io
#     locus_tag = ""
#     seq = ""
#     for ln in eachline(io)
#         if startswith(ln, '>')
#             locus_tag = split(split(ln, "locus_tag=")[2], ']'; limit = 2)[1]
#             if !haskey(seq_dict, locus_tag)
#                 seq_dict[locus_tag] = ""
#             end
#         else
#             seq = seq_dict[locus_tag]
#             seq_dict[locus_tag] = "$seq$ln"
#         end
#     end
# end

# CSV.write("data/databases/locus_tag_seq_ST398.csv", seq_dict; delim = '\t')

### make dic from chebi to smiles 
function read_json(file)
    open(file,"r") do f
        return JSON.parse(f)
    end
end
smiles_dic1 = read_json("data/databases/id-mapper1.json")
smiles_dic2 = read_json("data/databases/id-mapper2.json")
smiles_dic3 = read_json("data/databases/id-mapper3.json")

smiles_dic = merge(smiles_dic1,smiles_dic2,smiles_dic3)

chebi_2_smiles = Dict()
for (key,value) in smiles_dic
    for (x,y) in value
        if x == "SMILES"
            push!(chebi_2_smiles,key => y)
        end
    end
end

### make dict of chebi to inchi 
df = DataFrame(CSV.File("data/databases/locus_tag_seq_ST398.csv"))
g_id_sequence = Dict(Pair.(df.first, df.second))

chebi_2_inchi = DataFrame(CSV.File("data/databases/chebiId_inchi.tsv"))
chebi_2_inchi_dic = Dict(Pair.(chebi_2_inchi.CHEBI_ID, chebi_2_inchi.InChI))
chebi_2_inchi_dic = Dict("CHEBI:$x" => y for (x, y) in chebi_2_inchi_dic)

inchi_smiles = merge(chebi_2_smiles,chebi_2_inchi_dic)
### make dict of g_id => enyzme sequence 
df = DataFrame(
    reaction_id = String[],
    locus_tag = String[],
    amino_acid_sequence=String[],
    substrates=String[],
    products=String[],
)
problem_metabolites = DataFrame(mets = String[])
for (id, r) in model.reactions
    isnothing(r.gene_association_dnf) && continue
    subs = ""
    prods = ""
    for (x, y) in r.stoichiometry
        if !haskey(inchi_smiles,x)
            push!(problem_metabolites.mets,x)
            break 
        elseif y < 0
            subs = "$(inchi_smiles[x]);$subs"
        else
            prods = "$(inchi_smiles[x]);$prods"
        end
    end
    for iso in r.gene_association_dnf
        for g in iso
            push!(df, ("$(id)_f", g, g_id_sequence[g], subs, prods))
            push!(df, ("$(id)_r", g, g_id_sequence[g], prods, subs))
        end
    end
end

select!(df, [:amino_acid_sequence, :substrates, :products])

CSV.write("data/databases/reaction_sequence_subs_prods.csv", df; delim = ',')


CSV.write("problem_metabolites.csv", problem_metabolites)

### load reaction_isozymes

isozymes_df = DataFrame(CSV.File("data/model/reaction_isozymes.csv"))
select!(isozymes_df, :RHEA_ID, :Protein, :Stoichiometry, :Isozyme)

for (id, prot) in isozymes_df
    reaction_isozymes = Dict{String,Dict{String,Isozyme}}()

end

### change units of kcat if needed 

### add enzyme molar masses and change units if needed 

### fix total enzyme capacity

### run enzyme_constrained_flux_balance_analysis