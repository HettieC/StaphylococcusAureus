function enzyme_constraints!(model,reaction_isozymes)
    # add fake isozymes 
    oxphos_reactions = ["Ndh2", "Sdh", "Mqo", "Lqo", "cyt_aa3", "cyt_bd", "cyt_bo3","Ldh"]
    for rid in A.reactions(model)
        haskey(reaction_isozymes,rid) && continue
        if rid ∈ oxphos_reactions || !isnothing(tryparse(Int,rid))
            grrs = A.reaction_gene_association_dnf(model, rid)
            reaction_isozymes[rid] = Dict("isoyzme_$i" => Isozyme(
                gene_product_stoichiometry=Dict(g => 1 for g in grr), # assume subunit stoichiometry of 1 for all isozymes
                kcat_forward=65,
                kcat_reverse=65,
            )
            for (i,grr) in enumerate(model.reactions[rid].gene_association_dnf)
            )
        end
    end

    gene_product_molar_masses = get_gene_product_molar_mass([g for g in A.genes(model) if g != "g1"])
    gene_product_molar_masses["g1"] = sum(collect(values(gene_product_molar_masses)))/length(gene_product_molar_masses)

    # get membrane gids 
    membrane_gids = String[] 
    for (r,isos) in reaction_isozymes 
        if isnothing(tryparse(Int,r))
            append!(membrane_gids,[g for (x,y) in isos for (g,s) in y.gene_product_stoichiometry])
        end
    end
    df1 = DataFrame(CSV.File("data/databases/uniprot/s_aureus.csv"))
    df2 = DataFrame(CSV.File("data/databases/uniprot/idmapping_2025_05_19.csv"))
    rename!(df2,:Entry => :uniprot_accesion)
    df = innerjoin(df1,df2,on=:uniprot_accesion)
    subcellular_location = Dict{String,Vector{String}}()
    for ln in eachrow(df)        
        if ismissing(ln.Subcellular_location_CC)
            subcellular_location[ln.gene] = ["uncategorised"]
        else
            subs = string.(split(split(ln.Subcellular_location_CC,"SUBCELLULAR LOCATION: ")[2],"; "))
            subcellular_location[ln.gene] = rstrip.(first.(split.(subs," {")),'.')
        end
    end
    append!(membrane_gids, [x for (x,y) in subcellular_location if x ∈ A.genes(model) && any(z -> occursin("membrane",lowercase(z)),y) ])

    return gene_product_molar_masses, membrane_gids
end
export enzyme_constraints!
