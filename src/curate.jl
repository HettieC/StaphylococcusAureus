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
    end

    delete!(model.metabolites, "CHEBI:4167") # D-glucose
    delete!(model.metabolites, "CHEBI:61548") # D-glucose 6-phosphate
    delete!(model.metabolites, "CHEBI:91002") # sucrose 6(G)-phosphate

    # add generic transport gene
    model.genes["g1"] = CM.Gene(name="g1")

    # allow bidirectional H2O 
    model.reactions["EX_15377"].lower_bound = -1000
    model.reactions["EX_15377"].upper_bound = 1000

    #change directions to match what is found in biocyc - manual thermodynamics leaves much to be desired
    # biocyc = DataFrame(CSV.File(joinpath("data", "databases", "rhea", "biocyc_rxns.csv")))
    # @select!(biocyc, :rheaDir, :metacyc)
    # directions = String[]
    # for rid in A.reactions(model)
    #     isnothing(tryparse(Int,rid)) && continue
    #     qrt = RheaReactions.get_reaction_quartet(parse(Int, rid))
    #     df = @subset(biocyc, in.(:rheaDir, Ref(qrt)))
    #     isempty(df) && continue
    #     lb, ub = rhea_rxn_dir(df[1, 1], qrt)
    #     model.reactions[rid].lower_bound = lb
    #     model.reactions[rid].upper_bound = ub
    #     if lb != -1000 || ub != 1000
    #         push!(directions, rid)
    #     end
    # end

    #add atp maintenance reaction 
    model.reactions["ATPM"] = CM.Reaction(
        ;
        name="ATP maintenance",
        lower_bound=1.0,
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
            ]
        )
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
        )
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

            "CHEBI:46398" => -0.1,   #UTP
            "CHEBI:37565" => -0.059,   #GTP
             "CHEBI:37563" => -0.059,   #CTP 
             "CHEBI:57692" => -0.007,    #FAD  
            "CHEBI:61404" => -0.02,    #dATP
            "CHEBI:57287" => -4.42e-4, #CoA
             "CHEBI:37568" => -0.02,    #dTTP
             "CHEBI:61429" => -0.099,   #dGTP
             "CHEBI:61481" => -0.099,   #dCTP

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

            "CHEBI:57427" => -0.282,  #L-leucine
            "CHEBI:32682" => -0.111,  #L-arginine  
            "CHEBI:57762" => -0.207,  #L-valine  
            "CHEBI:60039" => -0.116,  #L-proline
            "CHEBI:35235" => -0.019,  #L-cysteine
            "CHEBI:57305" => -0.19,  #glycine
            "CHEBI:33384" => -0.19, #L-serine         
            "CHEBI:29991" => -0.261, #L-aspartate
            "CHEBI:57972" => -0.212, #L-alanine
            "CHEBI:58359" => -1.2,   #L-glutamine
            "CHEBI:29985" => -1.0,   #L-glutamate
            "CHEBI:32551" => -0.235, #L-lysine
            "CHEBI:58045" => -0.269, #L-isoleucine
            "CHEBI:57305" => -0.1,  #glycine
            "CHEBI:57926" => -0.179, #L-threonine
            "CHEBI:58095" => -0.137, #L-phenylalanine
            "CHEBI:58315" => -0.119, #L-tyrosine
            "CHEBI:57912" => -1.0,   #L-tryptophan
            "CHEBI:57595" => -0.073, #L-histidine
            "CHEBI:58199" => -0.1, #L-homocysteine
            "CHEBI:58048" => -0.1,  #L-asparagine
            "CHEBI:57844" => -0.084, #L-methionine

        ),
        objective_coefficient=1.0,
        notes=Dict("ref" => ["Diaz Calvo, S. epidermis, Metabolites 2022"]),
    )

    # add missing transporters 
    add_permease!(model, "CHEBI:32682", ["g1"], nothing)

    # make acsA bidirectional 
    model.reactions["23176"].lower_bound = -1000
    return model
end
