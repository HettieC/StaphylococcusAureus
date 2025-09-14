function add_electron_transport_chain!(model)

    model.metabolites["CHEBI:15378_p"] = CM.Metabolite(;
        name="H(+) periplasm",
        formula=Dict("H" => 1),
        charge=1,
        balance=0.0
    )

    model.reactions["Ndh2"] = CM.Reaction(;
        name="NADH dehydrogenase II",
        stoichiometry=Dict(
            "CHEBI:57945" => -1, #NADH 
            "CHEBI:44027" => -1, #menaquinone-8,
            "CHEBI:15378" => -1, #H+
            "CHEBI:57540" => 1, #NAD+ 
            "CHEBI:61684" => 1, #menaquinol-8 
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "NADH + menaquinone + H+ = NAD+ + menaquinol"
            ]
        ),
        gene_association_dnf = [
            ["SAPIG0228","SAPIG0518"]
        ]
    )

    model.reactions["Sdh"] = CM.Reaction(;
        name="Succinate dehydrogenase",
        stoichiometry=Dict(
            "CHEBI:44027" => -1, #menaquinone-8
            "CHEBI:30031" => -1, #succinate,
            "CHEBI:61684" => 1, #menaquinol-8 
            "CHEBI:29806" => 1, #fumarate 
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "menaquinone-8 + succinate = menaquinol-8 + fumarate"
            ]
        ),
        gene_association_dnf = [
            ["SAPIG1143"],
            #["SAPIG1144","SAPIG1145"]
        ]
    )

    model.reactions["Mqo"] = CM.Reaction(;
        name="Malate:quinone oxidoreductase",
        stoichiometry=Dict(
            "CHEBI:44027" => -1, #menaquinone-8 
            "CHEBI:15589" => -1, #malate,
            "CHEBI:61684" => 1, #menaquinol-8 
            "CHEBI:16452" => 1, #oxaloacetate 
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
            "CHEBI:44027" => -1, #menaquinone-8 
            "CHEBI:16651" => -1, #(S)-lactate
            "CHEBI:61684" => 1, #menaquinol-8 
            "CHEBI:15361" => 1, #pyruvate
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
            "CHEBI:57540" => -1, #NAD+ 
            "CHEBI:16651" => -1, #(S)-lactate
            "CHEBI:57945" => 1, #NADH 
            "CHEBI:15361" => 1, #pyruvate
            "CHEBI:15378" => 1, #H+
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
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -4, #H+
            "CHEBI:29033" => -4, #Fe(2+)
            "CHEBI:29034" => 4, #Fe(3+)
            "CHEBI:15377" => 2, #H2O
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
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -4, #H+
            "CHEBI:61684" => -2, #menaquinol-8
            "CHEBI:44027" => 2, #menaquinone-8
            "CHEBI:15377" => 2, #H2O
            "CHEBI:15378_p" => 4, #H+ periplasm
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "O2 + H+(c) + menaquinol = H2O + menaquinone + H+(p)"
            ]
        ),
        gene_association_dnf = [
            ["SAPIG1083","SAPIG1084"],
        ]
    )

    #cyt bo3
    model.reactions["cyt_bo3"] = CM.Reaction(;
        name="Cytochrome bo3 oxidase",
        stoichiometry=Dict(
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -8, #H+
            "CHEBI:61684" => -2, #menaquinol-8
            "CHEBI:44027" => 2, #menaquinone-8
            "CHEBI:15377" => 2, #H2O
            "CHEBI:15378_p" => 8, #H+ periplasm
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "O2 + H+(c) + menaquinol = H2O + menaquinone + H+(p)"
            ]
        ),
        gene_association_dnf = [
            ["g1"]
        ]
    )

    model
end
