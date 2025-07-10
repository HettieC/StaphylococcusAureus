function add_electron_transport_chain!(model)

    model.metabolites["15378_p"] = CM.Metabolite(;
        name="H(+) periplasm",
        formula=Dict("H" => 1),
        charge=1,
        balance=0.0
    )

    model.reactions["Ndh2"] = CM.Reaction(;
        name="NADH dehydrogenase II",
        stoichiometry=Dict(
            "57945" => -1, #NADH 
            "44027" => -1, #menaquinone-8,
            "15378" => -1, #H+
            "57540" => 1, #NAD+ 
            "61684" => 1, #menaquinol-8 
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "NADH + menaquinone-8 + H+ = NAD+ + menaquinol-8"
            ]
        ),
        gene_association_dnf = [
            ["SAPIG0228","SAPIG0518"]
        ]
    )

    model.reactions["Sdh"] = CM.Reaction(;
        name="Succinate dehydrogenase",
        stoichiometry=Dict(
            "57945" => -1, #NADH 
            "30031" => -1, #succinate,
            "15378" => -1, #H+
            "57540" => 1, #NAD+ 
            "29806" => 1, #fumarate 
            "15378_p" => 1, #H+ periplasm
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "NADH + succinate + H+ = NAD+ + fumarate + H+(p)"
            ]
        ),
        gene_association_dnf = [
            ["SAPIG1143"]
        ]
    )

    model.reactions["Mqo"] = CM.Reaction(;
        name="Malate:quinone oxidoreductase",
        stoichiometry=Dict(
            "44027" => -1, #menaquinone-8 
            "15589" => -1, #malate,
            "61684" => 1, #menaquinol-8 
            "16452" => 1, #oxaloacetate 
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "menaquinone-8 + malate = menaquinol-8 + oxaloacetate"
            ],
            "EC" => ["1.1.5.4"],
            "KEGG" => ["R01257","R00361"]
        ),
        gene_association_dnf = [
            ["SAPIG2418"],["SAPIG2655"]
        ]
    )    

    #Lqo 
    model.reactions["Lqo"] = CM.Reaction(;
        name="Lactate-quinone oxidoreductase",
        stoichiometry=Dict(
            "44027" => -1, #menaquinone-8 
            "16651" => -1, #(S)-lactate
            "61684" => 1, #menaquinol-8 
            "15361" => 1, #pyruvate
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "menaquinone-8 + (S)-lactate = menaquinol-8 + pyruvate"
            ],
        ),
        gene_association_dnf = [
            ["g1"]
        ]
    )

    model.reactions["Ldh"] = CM.Reaction(;
        name="Lactate dehydrogenase",
        stoichiometry=Dict(
            "57540" => -1, #NAD+ 
            "16651" => -1, #(S)-lactate
            "57945" => 1, #NADH 
            "15361" => 1, #pyruvate
            "15378" => 1, #H+
        ),
        annotations=Dict(
            "CM.Reaction" => [
                "NAD+ + (S)-lactate = NADH + pyruvate + H+"
            ],
            "EC" => ["1.1.1.27"],
            "KEGG" => ["R00703"]
        ),
        gene_association_dnf = [
            ["SAPIG0252","SAPIG2650"]
        ]
    )

    #cyt aa3 
    model.reactions["cyt_aa3"] = CM.Reaction(;
        name="Cytochrome aa3 oxidase",
        stoichiometry=Dict(
            "15379" => -1, #O2
            "15378" => -1, #H+
            "29033" => -1, #Fe(2+)
            "29034" => 1, #Fe(3+)
            "15377" => 1, #H2O
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "O2 + Fe(2+) + H+ = H2O + Fe(3+)"
            ]
        ),
        gene_association_dnf = [
            ["SAPIG1055","SAPIG1056","SAPIG1057","SAPIG1058"]
        ]
    )

    #cyt bd
    model.reactions["cyt_bd"] = CM.Reaction(;
        name="Cytochrome bd oxidase",
        stoichiometry=Dict(
            "15379" => -1, #O2
            "15378" => -4, #H+
            "61684" => -2, #menaquinol-8
            "44027" => 2, #menaquinone-8
            "15377" => 2, #H2O
            "15378_p" => 4, #H+ periplasm
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "O2 + H+(c) + menaquinol-8 = H2O + menaquinone-8 + H+(p)"
            ],
            "EC" => ["7.1.1.5"]
        ),
        gene_association_dnf = [
            ["SAPIG1083","SAPIG1084"],
            #["SAPIG1055","SAPIG1056","SAPIG1057"]
        ]
    )

    #cyt bo3
    model.reactions["cyt_bo3"] = CM.Reaction(;
        name="Cytochrome bo3 oxidase",
        stoichiometry=Dict(
            "15379" => -1, #O2
            "15378" => -8, #H+
            "61684" => -2, #menaquinol-8
            "44027" => 2, #menaquinone-8
            "15377" => 2, #H2O
            "15378_p" => 8, #H+ periplasm
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "O2 + H+(c) + menaquinol-8 = H2O + menaquinone-8 + H+(p)"
            ]
        ),
        gene_association_dnf = [
            ["g1"]
        ]
    )
    return model
end
