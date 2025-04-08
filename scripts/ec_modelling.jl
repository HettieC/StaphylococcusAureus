using COBREXA
using DataFrames, CSV

turnup_df = DataFrame(CSV.File("data/model/isozymes/turnup_output.csv")) #forward and reverse reactions have the same kcats?

rxn_seq_subs_prods = DataFrame(CSV.File("data/model/isozymes/reaction_sequence_subs_prods.csv"))

df = DataFrame(CSV.File("data/model/isozymes/locus_tag_seq_ST398.csv"))
g_id_sequence = Dict(Pair.(df.first, df.second))

turnup_df = DataFrames.rename!(turnup_df,
    "kcat [s^(-1)]" => "kcat",
)
kcat_df = insertcols(rxn_seq_subs_prods, 5, :kcat => turnup_df.kcat)
kcat_dict = Dict{String,Dict{String,Float64}}()
for row in eachrow(kcat_df)
    if !ismissing(row.kcat)
        if !haskey(kcat_dict, row.reaction_id)
            # put into 1/h instead of 1/s by multiplying by 3.6
            kcat_dict[row.reaction_id] = Dict(row.locus_tag => row.kcat * 3.6)
        else
            kcat_dict[row.reaction_id][row.locus_tag] = row.kcat * 3.6
        end
    end
end
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

# calculate molar masses of gene products

function get_gene_product_molar_mass(gids)
    AA_mass = Dict(
        'A' => 89,
        'R' => 174,
        'N' => 132,
        'D' => 133,
        'B' => 133,
        'C' => 121,
        'Q' => 146,
        'E' => 147,
        'Z' => 147,
        'G' => 75,
        'H' => 155,
        'I' => 131,
        'L' => 131,
        'K' => 146,
        'M' => 149,
        'F' => 165,
        'P' => 115,
        'S' => 105,
        'T' => 119,
        'W' => 204,
        'Y' => 181,
        'V' => 117,
    )
    seqs = Dict{String,String}()
    open("data/sequence.txt","r") do io
        prev = ""
        for ln in eachline(io)
            if startswith(ln,'>')
                seqs[ln[2:end]] = ""
                prev = ln[2:end]
            else 
                seqs[prev] = seqs[prev]*ln 
            end
        end
    end

    gene_product_molar_mass = Dict{String,Float64}()
    for gid in gids 
        gene_product_molar_mass[gid] = sum([AA_mass[aa] for aa in seqs[gid]]) / 1000 # convert to kDa
    end

    println(110*length(seqs["SAPIG1968"]))
    return gene_product_molar_mass     
    
end

molar_masses = get_gene_product_molar_mass(A.genes(model))

# convert units from g/mol to kg/mol

# for (met,mass) in molar_masses
#     molar_masses[met] = molar_masses[met] / 1000
# end


total_enzyme_capacity = 50.0 # mg of enzyme/gDW

ec_sol = ec_solution = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses = molar_masses,
    capacity = total_enzyme_capacity,
    optimizer = HiGHS.Optimizer,
)

open("data/fluxes.json", "w") do io
    JSON.print(io, Dict(string(x) => y for (x, y) in sol.fluxes))
end

# total amount of required gene product mass
ec_solution.gene_product_capacity

# amount of gene product material required for the system to run
gene_product_amounts = ec_solution.gene_product_amounts

# FBA
sol = parsimonious_flux_balance_analysis(model; optimizer=HiGHS.Optimizer)
