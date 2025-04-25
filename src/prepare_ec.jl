using DataFrames, CSV
using CairoMakie
using Statistics
using Serialization

turnup_df = DataFrame(CSV.File("data/model/isozymes/turnup_output.csv"))

rxn_seq_subs_prods = DataFrame(CSV.File("data/model/isozymes/reaction_sequence_subs_prods.csv"))

df = DataFrame(CSV.File("data/model/isozymes/locus_tag_seq_ST398.csv"))
g_id_sequence = Dict(Pair.(df.first, df.second))

turnup_df = DataFrames.rename!(turnup_df,
    "kcat [s^(-1)]" => "kcat",
)
kcat_df = insertcols(rxn_seq_subs_prods, 5, :kcat => turnup_df.kcat)
# kcat mean at 22.82
kcat_dict = Dict{String,Dict{String,Float64}}()
for row in eachrow(kcat_df)
    if !ismissing(row.kcat)
        if !haskey(kcat_dict, row.reaction_id)
            # put into 1/h instead of 1/s by multiplying by 3600
            kcat_dict[row.reaction_id] = Dict(row.locus_tag => row.kcat * 3.6)
        else
            kcat_dict[row.reaction_id][row.locus_tag] = row.kcat * 3.6
        end
    end
end
avg_kcat = mean(vcat([collect(values(y)) for (x,y) in kcat_dict]...))

isozymes_stoich_df = DataFrame(CSV.File("data/model/isozymes/reaction_isozymes.csv")) 
select!(isozymes_stoich_df,:RHEA_ID, :Protein, :Stoichiometry, :Isozyme)
isozyme_stoich = Dict(Pair.(isozymes_stoich_df.Protein, isozymes_stoich_df.Stoichiometry))

reaction_isozymes = Dict{String,Dict{String,Isozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for (id, rxn) in model.reactions
    grrs = rxn.gene_association_dnf
    isnothing(grrs) && continue # skip if no grr available
    haskey(kcat_dict, "$(id)_f") || continue # skip if no kcat data available
    haskey(kcat_dict, "$(id)_r") || continue
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, id, Dict{String,Isozyme}())
        #println(id,"  ",[isozyme_stoich[g] for g in grr])
        d["isozyme_"*string(i)] = Isozyme( # each isozyme gets a unique name
            gene_product_stoichiometry = Dict(g => isozyme_stoich[g] for g in grr), 
            kcat_forward = maximum([kcat_dict["$(id)_f"][g] for g in grr]),
            kcat_reverse = maximum([kcat_dict["$(id)_r"][g] for g in grr]),
        )
    end
end


# add oxphos fake isozymes 
oxphos_reactions = ["Ndh2", "Sdh", "Mqo", "Lqo", "cyt_aa3", "cyt_bd", "cyt_bo3"]
for rid in oxphos_reactions
    grrs = A.reaction_gene_association_dnf(model, rid)
    println(grrs)
    reaction_isozymes[rid] = Dict("isoyzme_1" => Isozyme(
        gene_product_stoichiometry = Dict("g1" => 1), # assume subunit stoichiometry of 1 for all isozymes
        kcat_forward = avg_kcat,
        kcat_reverse = avg_kcat,
        )
    )
end

serialize("reaction_isozymes.serialized", reaction_isozymes)