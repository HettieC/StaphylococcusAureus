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
            "POLYMER:9537" => -1, #a menaquinone,
            "CHEBI:15378" => -1, #H+
            "CHEBI:57540" => 1, #NAD+ 
            "POLYMER:9539" => 1, #a menaquinol 
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "NADH + menaquinone + H+ = NAD+ + menaquinol"
            ]
        )
    )

    model.reactions["Sdh"] = CM.Reaction(;
        name="Succinate dehydrogenase",
        stoichiometry=Dict(
            "CHEBI:57945" => -1, #NADH 
            "CHEBI:30031" => -1, #succinate,
            "CHEBI:15378" => -1, #H+
            "CHEBI:57540" => 1, #NAD+ 
            "CHEBI:29806" => 1, #fumarate 
            "CHEBI:15378_p" => 1, #H+ periplasm
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "NADH + succinate + H+ = NAD+ + fumarate + H+(p)"
            ]
        )
    )

    model.reactions["Mqo"] = CM.Reaction(;
        name="Malate:quinone oxidoreductase",
        stoichiometry=Dict(
            "CHEBI:57945" => -1, #NADH 
            "CHEBI:15589" => -1, #malate,
            "CHEBI:57540" => 1, #NAD+ 
            "CHEBI:16452" => 1, #oxaloacetate 
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "NADH + malate = NAD+ + oxaloacetate"
            ]
        )
    )

    #Lqo 
    model.reactions["Lqo"] = CM.Reaction(;
        name="Lactate-quinone oxidoreductase",
        stoichiometry=Dict(
            "POLYMER:9537" => -1, #a menaquinone 
            "CHEBI:16651" => -1, #(S)-lactate
            "POLYMER:9539" => 1, #a menaquinol 
            "CHEBI:15361" => 1, #pyruvate
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "a menaquinone + (S)-lactate = a menaquinol + pyruvate"
            ],
        )
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
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "NAD+ + (S)-lactate = NADH + pyruvate + H+"
            ],
        )
    )

    #cyt aa3 
    model.reactions["cyt_aa3"] = CM.Reaction(;
        name="Cytochrome aa3 oxidase",
        stoichiometry=Dict(
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -1, #H+
            "CHEBI:29033" => -1, #Fe(2+)
            "CHEBI:29034" => 1, #Fe(3+)
            "CHEBI:15377" => 1, #H2O
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "O2 + Fe(2+) + H+ = H2O + Fe(3+)"
            ]
        )
    )

    #cyt bd
    model.reactions["cyt_bd"] = CM.Reaction(;
        name="Cytochrome bd oxidase",
        stoichiometry=Dict(
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -1, #H+
            "POLYMER:9539" => -1, #a menaquinol
            "POLYMER:9537" => 1, #a menaquinone
            "CHEBI:15377" => 1, #H2O
            "CHEBI:15378_p" => 1, #H+ periplasm
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "O2 + H+(c) + menaquinol = H2O + menaquinone + H+(p)"
            ]
        )
    )

    #cyt bo3
    model.reactions["cyt_bo3"] = CM.Reaction(;
        name="Cytochrome bo3 oxidase",
        stoichiometry=Dict(
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -1, #H+
            "POLYMER:9539" => -1, #a menaquinol
            "POLYMER:9537" => 1, #a menaquinone
            "CHEBI:15377" => 1, #H2O
            "CHEBI:15378_p" => 1, #H+ periplasm
        ),
        lower_bound = 0.0,
        annotations=Dict(
            "CM.Reaction" => [
                "O2 + H+(c) + menaquinol = H2O + menaquinone + H+(p)"
            ]
        )
    )

    model
end
