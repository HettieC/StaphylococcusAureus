using StaphylococcusAureus
import AbstractFBCModels as A
import ConstraintTrees as C
using COBREXA
using DataFrames, CSV
using CairoMakie
using Statistics
using HiGHS


model = build_model()

reaction_isozymes = get_reaction_isozymes(model)

avg_kcat = mean(vcat([b.kcat_forward for (x,y) in reaction_isozymes for (a,b) in y]...))

# add oxphos fake isozymes 
oxphos_reactions = ["Ndh2", "Sdh", "Mqo", "Lqo", "cyt_aa3", "cyt_bd", "cyt_bo3"]
for rid in oxphos_reactions
    grrs = A.reaction_gene_association_dnf(model, rid)
    reaction_isozymes[rid] = Dict("isoyzme_1" => Isozyme(
        gene_product_stoichiometry = Dict("g1" => 1), # assume subunit stoichiometry of 1 for all isozymes
        kcat_forward = avg_kcat,
        kcat_reverse = avg_kcat,
        )
    )
end


gene_product_molar_masses = get_gene_product_molar_mass([g for g in A.genes(model) if g!="g1"])
gene_product_molar_masses["g1"] = mean(collect(values(gene_product_molar_masses)))


total_capacity = 500.0 #


ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses = gene_product_molar_masses,
    capacity = total_capacity,
    optimizer = HiGHS.Optimizer,
)


for (r,rxn) in model.reactions 
    haskey(reaction_isozymes,r) && continue 
    isnothing(tryparse(Int64,r)) && continue
    for ec in rxn.annotations["EC"]
        df = find_ec_info(ec)
        isempty(df) && continue 
        println(rxn.annotations["KEGG"],"  ", r, "  ", rxn.annotations["EC"])
        println(df)
    end
end

