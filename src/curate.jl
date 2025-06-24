rhea_rxn_dir(rxn, qrt) = begin
    idx = first(indexin([rxn], qrt))
    isnothing(idx) && error("Reaction not found...")
    idx == 1 && return (-1000, 1000)
    idx == 2 && return (0, 1000)
    idx == 3 && return (-1000, 0)
    idx == 4 && return (-1000, 1000)
end

function curate!(model)

    # modify rhea reactions to use beta-D isomer instead of D-glucose 
    # modify rhea reactions to use general sucrose 6(F)-phosphate instead of 6G
    # modify rhea reactions to use beta-D isomer instead of D-fructose 
    # modify rhea reactions to use menaquinol/menaquinone instead of ubiquinol/ubiquinone
    # modify rhea reactions to use (S)-1-pyrroline-5-carboxylate instead of 1-pyrroline-5-carboxylate 


    # #remove general quinone/quinol/ubiquinone/ubiquinone reactions 
    for (r,rxn) in model.reactions 
        if haskey(rxn.stoichiometry,"CHEBI:132124") || haskey(rxn.stoichiometry,"CHEBI:24646") || haskey(rxn.stoichiometry,"POLYMER:9566") || haskey(rxn.stoichiometry,"POLYMER:9565") || haskey(rxn.stoichiometry,"POLYMER:9563")
            delete!(model.reactions,r)
        end
    end

    for (r, rxn) in model.reactions
        if haskey(rxn.stoichiometry, "CHEBI:4167") # D-glucose -> beta-D-glucose
            rxn.stoichiometry["CHEBI:15903"] = rxn.stoichiometry["CHEBI:4167"]
            delete!(rxn.stoichiometry, "CHEBI:4167")
        end
        if haskey(rxn.stoichiometry, "CHEBI:61548") # D-glucose 6-phosphate -> β-D-glucose 6-phosphate 
            rxn.stoichiometry["CHEBI:58247"] = rxn.stoichiometry["CHEBI:61548"]
            delete!(rxn.stoichiometry, "CHEBI:61548")
        end
        if haskey(rxn.stoichiometry, "CHEBI:91002") # sucrose 6(G)-phosphate -> sucrose 6(F)-phosphate 
            rxn.stoichiometry["CHEBI:57723"] = rxn.stoichiometry["CHEBI:91002"]
            delete!(rxn.stoichiometry, "CHEBI:91002")
        end
        if haskey(rxn.stoichiometry, "CHEBI:37721") # D-fructose -> β-D-fructose
            rxn.stoichiometry["CHEBI:28645"] = rxn.stoichiometry["CHEBI:37721"]
            delete!(rxn.stoichiometry, "CHEBI:37721")
        end
        if haskey(rxn.stoichiometry, "CHEBI:61527") # D-fructose -> β-D-fructose
            rxn.stoichiometry["CHEBI:57634"] = rxn.stoichiometry["CHEBI:61527"]
            delete!(rxn.stoichiometry, "CHEBI:61527")
        end
        if haskey(rxn.stoichiometry, "CHEBI:16389") # a ubiquinone -> a menaquinone
            rxn.stoichiometry["POLYMER:9537"] = rxn.stoichiometry["CHEBI:16389"]
            delete!(rxn.stoichiometry, "CHEBI:16389")
        end
        if haskey(rxn.stoichiometry, "CHEBI:17976") # a ubiquinol -> a menaquinol
            rxn.stoichiometry["POLYMER:9539"] = rxn.stoichiometry["CHEBI:17976"]
            delete!(rxn.stoichiometry, "CHEBI:17976")
        end
        if haskey(rxn.stoichiometry, "CHEBI:15893") # 1-pyrroline-5-carboxylate  -> (S)-1-pyrroline-5-carboxylate 
            rxn.stoichiometry["CHEBI:17388"] = rxn.stoichiometry["CHEBI:15893"]
            delete!(rxn.stoichiometry, "CHEBI:15893")
        end
    end


    delete!(model.metabolites, "CHEBI:4167") # D-glucose
    delete!(model.metabolites, "CHEBI:61548") # D-glucose 6-phosphate
    delete!(model.metabolites, "CHEBI:91002") # sucrose 6(G)-phosphate

    # add generic transport gene
    model.genes["g1"] = CM.Gene(name="g1")

    # allow bidirectional H2O 
    model.reactions["EX_15377"].lower_bound = -1000
    model.reactions["EX_15377"].upper_bound = 1000

    #change directions to match what is found in biocyc 
    biocyc = DataFrame(CSV.File(joinpath("data", "databases", "rhea", "biocyc_rxns.csv")))
    bidirectional = string.(JSON.parsefile("data/model/reactions/bidirectional.json"))
    manual_directions = DataFrame(CSV.File("data/model/reactions/unidirectional_reactions.csv"))
    @select!(biocyc, :rheaDir, :metacyc)
    for rid in A.reactions(model)
        rid ∈ bidirectional && continue
        isnothing(tryparse(Int,rid)) && continue
        rid ∈ string.(manual_directions.RHEA_ID) && continue # ignore if direction manually specified
        qrt = RheaReactions.get_reaction_quartet(parse(Int, rid))
        df = @subset(biocyc, in.(:rheaDir, Ref(qrt)))
        isempty(df) && continue
        lb, ub = rhea_rxn_dir(df[1, 1], qrt)
        model.reactions[rid].lower_bound = lb
        model.reactions[rid].upper_bound = ub
    end

    #add atp maintenance reaction 
    model.reactions["ATPM"] = CM.Reaction(
        ;
        name="ATP maintenance",
        lower_bound=5.0,
        stoichiometry=Dict(
            "CHEBI:30616" => -1, #atp
            "CHEBI:15377" => -1, #h2o
            "CHEBI:43474" => 1, #phosphate
            "CHEBI:15378" => 1, #h+
            "CHEBI:456216" => 1, #adp
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
            "CHEBI:456216" => -1, #adp
            "CHEBI:15378_p" => -4, #h+ periplasm
            "CHEBI:43474" => -1, #phosphate
            "CHEBI:30616" => 1, #atp
            "CHEBI:15377" => 1, #h2o
            "CHEBI:15378" => 3, #h+
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
            "CHEBI:30616" => -50, #atp
            "CHEBI:15377" => -50, #h2o
            "CHEBI:43474" => 50, #phosphate
            "CHEBI:15378" => 50, #h+
            "CHEBI:456216" => 50, #adp

            "CHEBI:46398" => -0.506,   #UTP
            "CHEBI:37565" => -0.496,   #GTP
            "CHEBI:37563" => -0.496,   #CTP 
            "CHEBI:57692" => -0.00067,    #FAD  
            "CHEBI:61404" => -0.676,    #dATP
            "CHEBI:57287" => -4.42e-4, #CoA
            "CHEBI:37568" => -0.676,    #dTTP
            "CHEBI:61429" => -0.33,   #dGTP
            "CHEBI:61481" => -0.33,   #dCTP

            "CHEBI:57288" => -0.00033, #acetyl-coa

            "CHEBI:57783" => 2e-4, #NADPH 
            "CHEBI:57945" => 2e-4, #NADH
            "CHEBI:58349" => -2e-4, #NADP(+)
            "CHEBI:57540" => -2e-4, #NAD(+)

            "CHEBI:30807" => -0.1,    #tetradecanoate
            "CHEBI:25646" => -0.1,    #octanoate
            "CHEBI:7896" => -0.1,     #hexadecanoate
            "CHEBI:18262" => -0.1,    #dodecanoate
            "CHEBI:27689" => -0.1,    #decanoate
            "CHEBI:25629" => -0.1,    #octadecanoate
            # "POLYMER:9537" => -0.0765, #menaquinone 
            # "POLYMER:9539" => -0.0765, #menaquinol

            #protein
            "CHEBI:57427" => -0.699,  #L-leucine
            "CHEBI:32682" => -0.268,  #L-arginine  
            "CHEBI:57762" => -0.514,  #L-valine  
            "CHEBI:60039" => -0.246,  #L-proline
            "CHEBI:35235" => -0.048,  #L-cysteine
            "CHEBI:57305" => -0.19,  #glycine
            "CHEBI:33384" => -0.469, #L-serine         
            "CHEBI:29991" => -0.446, #L-aspartate
            "CHEBI:57972" => -0.492, #L-alanine
            "CHEBI:58359" => -0.318,  #L-glutamine
            "CHEBI:29985" => -0.497,   #L-glutamate
            "CHEBI:32551" => -0.576, #L-lysine
            "CHEBI:58045" => -0.658, #L-isoleucine
            "CHEBI:57305" => -0.462,  #glycine
            "CHEBI:57926" => -0.443, #L-threonine
            "CHEBI:58095" => -0.246, #L-phenylalanine
            "CHEBI:58315" => -0.298, #L-tyrosine
            "CHEBI:57912" => -0.057,   #L-tryptophan
            "CHEBI:57595" => -0.178, #L-histidine
            "CHEBI:58199" => -0.1, #L-homocysteine
            "CHEBI:58048" => -0.433,  #L-asparagine
            "CHEBI:57844" => -0.202, #L-methionine
        ),
        objective_coefficient=1.0,
        notes=Dict("ref" => ["Diaz Calvo, S. epidermis, Metabolites 2022"]),
    )

    # add missing transporters 
    add_permease!(model, "CHEBI:32682", ["g1"], nothing)
    return model
end

