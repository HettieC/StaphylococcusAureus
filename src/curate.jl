rhea_rxn_dir(rid, consensus) = begin
    idx = consensus - parse(Int, rid)
    idx == 0 && return (-1000, 1000)
    idx == 1 && return (0, 1000)
    idx == 2 && return (-1000, 0)
    idx == 3 && return (-1000, 1000)
end

function curate!(model)

    # modify rhea reactions to use beta-D isomer instead of D-glucose 
    # modify rhea reactions to use general sucrose 6(F)-phosphate instead of 6G
    # modify rhea reactions to use beta-D isomer instead of D-fructose 
    # modify rhea reactions to use menaquinol/menaquinone instead of ubiquinol/ubiquinone
    # modify rhea reactions to use (S)-1-pyrroline-5-carboxylate instead of 1-pyrroline-5-carboxylate 


    # #remove general quinone/quinol/ubiquinone/ubiquinone reactions 
    for (r,rxn) in model.reactions 
        if haskey(rxn.stoichiometry,"132124") || haskey(rxn.stoichiometry,"24646") || haskey(rxn.stoichiometry,"POLYMER:9566") || haskey(rxn.stoichiometry,"POLYMER:9565") || haskey(rxn.stoichiometry,"POLYMER:9563")
            delete!(model.reactions,r)
        end
    end

    for (r, rxn) in model.reactions
        if haskey(rxn.stoichiometry, "4167") # D-glucose -> beta-D-glucose
            rxn.stoichiometry["15903"] = rxn.stoichiometry["4167"]
            delete!(rxn.stoichiometry, "4167")
        end
        if haskey(rxn.stoichiometry, "61548") # D-glucose 6-phosphate -> β-D-glucose 6-phosphate 
            rxn.stoichiometry["58247"] = rxn.stoichiometry["61548"]
            delete!(rxn.stoichiometry, "61548")
        end
        if haskey(rxn.stoichiometry, "91002") # sucrose 6(G)-phosphate -> sucrose 6(F)-phosphate 
            rxn.stoichiometry["57723"] = rxn.stoichiometry["91002"]
            delete!(rxn.stoichiometry, "91002")
        end
        if haskey(rxn.stoichiometry, "37721") # D-fructose -> β-D-fructose
            rxn.stoichiometry["28645"] = rxn.stoichiometry["37721"]
            delete!(rxn.stoichiometry, "37721")
        end
        if haskey(rxn.stoichiometry, "61527") # D-fructose -> β-D-fructose
            rxn.stoichiometry["57634"] = rxn.stoichiometry["61527"]
            delete!(rxn.stoichiometry, "61527")
        end
        if haskey(rxn.stoichiometry, "16389") # a ubiquinone -> a menaquinone
            rxn.stoichiometry["POLYMER:9537"] = rxn.stoichiometry["16389"]
            delete!(rxn.stoichiometry, "16389")
        end
        if haskey(rxn.stoichiometry, "17976") # a ubiquinol -> a menaquinol
            rxn.stoichiometry["POLYMER:9539"] = rxn.stoichiometry["17976"]
            delete!(rxn.stoichiometry, "17976")
        end
        if haskey(rxn.stoichiometry, "15893") # 1-pyrroline-5-carboxylate  -> (S)-1-pyrroline-5-carboxylate 
            rxn.stoichiometry["17388"] = rxn.stoichiometry["15893"]
            delete!(rxn.stoichiometry, "15893")
        end
    end


    delete!(model.metabolites, "4167") # D-glucose
    delete!(model.metabolites, "61548") # D-glucose 6-phosphate
    delete!(model.metabolites, "91002") # sucrose 6(G)-phosphate

    # add generic transport gene
    model.genes["g1"] = CM.Gene(name="g1")

    #change directions to match what is found in biocyc 
    # bidirectional = string.(JSON.parsefile("data/model/reactions/bidirectional.json"))
    # manual_directions = DataFrame(CSV.File("data/model/reactions/unidirectional_reactions.csv"))
    # metacyc = DataFrame(
    #     CSV.File(
    #         joinpath(pkgdir(@__MODULE__), "data", "databases", "rhea", "biocyc_rxns.csv"),
    #         drop = [3],
    #     ),
    # )
    # @rename!(metacyc, :metacyc = :rheaDir)

    # ecocyc = DataFrame(
    #     CSV.File(
    #         joinpath(pkgdir(@__MODULE__), "data", "databases", "rhea", "ecocyc_rxns.csv"),
    #         drop = [3],
    #     ),
    # )
    # @rename!(ecocyc, :ecocyc = :rheaDir)

    # j = outerjoin(metacyc, ecocyc, on = :rhea)
    # @rtransform!(
    #     j,
    #     :consensus = ismissing(:ecocyc) ? :metacyc : :ecocyc,
    #     :rhea = string(:rhea)
    # )
    # for rid in A.reactions(model)
    #     rid ∈ bidirectional && continue
    #     isnothing(tryparse(Int,rid)) && continue
    #     rid ∈ string.(manual_directions.RHEA_ID) && continue # ignore if direction manually specified
    #     #TODO! speed up get_quartet by using rhea database locally!
    #     qrt = get_quartet(parse(Int, rid))
    #     master_rid = collect(keys(qrt))[1]
    #     df = @subset(j, in.(:rhea, Ref(qrt[master_rid])))
    #     isempty(df) && continue
    #     lb, ub = rhea_rxn_dir(master_rid, df.consensus[1])
    #     model.reactions[rid].lower_bound = lb
    #     model.reactions[rid].upper_bound = ub
    # end

    #add atp maintenance reaction 
    model.reactions["ATPM"] = CM.Reaction(
        ;
        name="ATP maintenance",
        lower_bound=5.0,
        stoichiometry=Dict(
            "30616" => -1, #atp
            "15377" => -1, #h2o
            "43474" => 1, #phosphate
            "15378" => 1, #h+
            "456216" => 1, #adp
        ),
        annotations = Dict(
            "CM.Reaction" => [
                "ATP + H2O = ADP + phosphate + H+"
            ],
            "EC" => ["3.6.1.5 3.6.1.8"],
            "KEGG" => ["R00086"]
        ),
        gene_association_dnf = [["SAPIG2145"],["SAPIG2147"]],
    )

    #add atp synthase reaction 
    model.reactions["ATPS"] = CM.Reaction(
        ;
        name="ATP synthase",
        lower_bound=0.0,
        stoichiometry=Dict(
            "456216" => -1, #adp
            "15378_p" => -4, #h+ periplasm
            "43474" => -1, #phosphate
            "30616" => 1, #atp
            "15377" => 1, #h2o
            "15378" => 3, #h+
        ),
        annotations = Dict(
            "CM.Reaction" => [
                "ADP + phosphate + 4 H+ (periplasm) = ATP + H2O + 3 H+"
            ]
        ),
        gene_association_dnf = [["SAPIG2145"], ["SAPIG2147"]]
    )

    # add a biomass reaction
    model.reactions["biomass"] = CM.Reaction(
        ;
        name="Biomass based on Staph epidermis RP62A",
        lower_bound=0.0,
        upper_bound=1000.0,
        stoichiometry=Dict(
            "30616" => -50, #atp
            "15377" => -50, #h2o
            "43474" => 50, #phosphate
            "15378" => 50, #h+
            "456216" => 50, #adp

            "46398" => -0.506,   #UTP
            "37565" => -0.496,   #GTP
            "37563" => -0.496,   #CTP 
            "57692" => -0.00067,    #FAD  
            "61404" => -0.676,    #dATP
            "57287" => -4.42e-4, #CoA
            "37568" => -0.676,    #dTTP
            "61429" => -0.33,   #dGTP
            "61481" => -0.33,   #dCTP

            "57288" => -0.00033, #acetyl-coa

            "57783" => 2e-4, #NADPH 
            "57945" => 2e-4, #NADH
            "58349" => -2e-4, #NADP(+)
            "57540" => -2e-4, #NAD(+)

            "30807" => -0.1,    #tetradecanoate
            "25646" => -0.1,    #octanoate
            "7896" => -0.1,     #hexadecanoate
            "18262" => -0.1,    #dodecanoate
            "27689" => -0.1,    #decanoate
            "25629" => -0.1,    #octadecanoate
            # "POLYMER:9537" => -0.0765, #menaquinone 
            # "POLYMER:9539" => -0.0765, #menaquinol

            #protein
            "57427" => -0.699,  #L-leucine
            "32682" => -0.268,  #L-arginine  
            "57762" => -0.514,  #L-valine  
            "60039" => -0.246,  #L-proline
            "35235" => -0.048,  #L-cysteine
            "57305" => -0.19,  #glycine
            "33384" => -0.469, #L-serine         
            "29991" => -0.446, #L-aspartate
            "57972" => -0.492, #L-alanine
            "58359" => -0.318,  #L-glutamine
            "29985" => -0.497,   #L-glutamate
            "32551" => -0.576, #L-lysine
            "58045" => -0.658, #L-isoleucine
            "57305" => -0.462,  #glycine
            "57926" => -0.443, #L-threonine
            "58095" => -0.246, #L-phenylalanine
            "58315" => -0.298, #L-tyrosine
            "57912" => -0.057,   #L-tryptophan
            "57595" => -0.178, #L-histidine
            "58199" => -0.1, #L-homocysteine
            "58048" => -0.433,  #L-asparagine
            "57844" => -0.202, #L-methionine
        ),
        objective_coefficient=1.0,
        notes=Dict("ref" => ["Diaz Calvo, S. epidermis, Metabolites 2022"]),
    )

    # add missing transporters 
    add_permease!(model, "32682", ["g1"], nothing)
    return model
end

