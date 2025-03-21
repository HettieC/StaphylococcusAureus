using COBREXA
using DataFrames, CSV

turnup_df = DataFrame(CSV.File("data/databases/turnup_output.csv"))

rxn_seq_subs_prods = DataFrame(CSV.File("data/databases/reaction_sequence_subs_prods.csv"))

df = DataFrame(CSV.File("data/databases/locus_tag_seq_ST398.csv"))
g_id_sequence = Dict(Pair.(df.first, df.second))

turnup_df = DataFrames.rename!(turnup_df,
    "kcat [s^(-1)]" => "kcat",
)
kcat_df = insertcols(rxn_seq_subs_prods, 5, :kcat => turnup_df.kcat)
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
isozymes_stoich_df = DataFrame(CSV.File("data/model/reaction_isozymes.csv")) 
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
            gene_product_stoichiometry = Dict(g => isozyme_stoich[g] for g in grr), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = maximum([kcat_dict["$(id)_f"][g] for g in grr]),
            kcat_reverse = maximum([kcat_dict["$(id)_r"][g] for g in grr]),
        )
    end
end

molar_masses = Dict{String,Float64}() # g/mol
begin
    molar_masses["61404"] = 487.1499 # dATP
    molar_masses["61429"] = 503.1493 # dGTP
    molar_masses["61481"] = 463.1252 # dCTP
    molar_masses["37568"] = 478.1365 # dTTP

    molar_masses["30616"] = 503.14946 # ATP
    molar_masses["37565"] = 519.14886 # GTP
    molar_masses["37563"] = 479.12468 # CTP
    molar_masses["46398"] = 480.1094 # UTP

    molar_masses["32551"] = 147.19558 # # lysine
    molar_masses["58045"] = 131.1729 # isoleucine
    molar_masses["57427"] = 131.1729 # leucine
    molar_masses["57844"] = 149.2124 # methionine
    molar_masses["58095"] = 165.1891 # phenylalanine
    molar_masses["57926"] = 119.1197  # threonine
    molar_masses["57912"] = 204.2262 # tryptophan
    molar_masses["57762"] = 117.1469 # valine
    molar_masses["32682"] = 175.20906 # arginine
    molar_masses["57595"] = 155.1552 # histidine
    molar_masses["57972"] = 89.0935 # alanine
    molar_masses["58048"] = 132.1184 # asparagine
    molar_masses["29991"] = 132.09478 # aspartate
    molar_masses["35235"] = 121.159 # cysteine
    molar_masses["29985"] = 147.1299 # glutamate
    molar_masses["58359"] = 146.1451 # glutamine
    molar_masses["57305"] = 75.0669 # glycine
    molar_masses["60039"] = 115.131 # proline
    molar_masses["33384"] = 105.093 # serine
    molar_masses["58315"] = 181.1894 # tyrosine

    molar_masses["glycogen"] = 162.1406 # C6H10O5

    molar_masses["30807"] = 227.364 # tetradecanoic acid
    molar_masses["7896"] = 255.4161 # hexadecanoic acid
    molar_masses["25629"] = 283.47 # octadecanoic acid

    molar_masses["peptidoglycan"] = 1916.20990

    molar_masses["kdo_lps"] = 2232.67080

    # soluble pool
    molar_masses["60530"] = 836.838 # Fe(II)-heme o
    molar_masses["57692"] = 782.5259 # FAD
    molar_masses["57705"] = 605.3378 # UDP-N-acetyl-alpha-D-glucosamine
    molar_masses["57540"] = 662.4172 # NAD(+)
    molar_masses["58885"] = 564.2859 # UDP-alpha-D-glucose
    molar_masses["57287"] = 763.502 # CoA
    molar_masses["57925"] = 306.31 # glutathione
    molar_masses["57945"] = 663.4251 # NADH
    molar_masses["58223"] = 401.1374 # UDP
    molar_masses["29985"] = 146.12136 # L-glutamate
    molar_masses["32966"] = 336.08392 # beta-D-fructose 1,6-bisphosphate
    molar_masses["30616"] = 503.14946 # ATP
    molar_masses["57783"] = 741.3891 # NADPH
    molar_masses["57986"] = 375.356 # riboflavin
    molar_masses["597326"] = 245.126 # pyridoxal 5'-phosphate
    molar_masses["62501"] = 439.3816 # folate
    molar_masses["58297"] = 610.615 # glutathione disulfide
    molar_masses["58210"] = 453.325 # FMN
    molar_masses["58349"] = 740.3812 # NADP(+)
end

total_enzyme_capacity = 50.0 # mg of enzyme/gDW

ec_solution = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses = molar_masses,
    capacity = total_enzyme_capacity,
    optimizer = HiGHS.Optimizer,
)