function add_electron_transport_chain!(model)

    gs = String[] 
    ms = String[] 

    # Ndh2
    rhea_ids = [
        "29103",#MQ8
        "35559",#MQ7
    #Sdh
        "71663",#general proton-pumping
        "27834",#general non-H+ pumping
        "29215",#MQ8 non-H+ pumping
        "33863",#Mq7 non-H+ pumping
    #Mqo
        "30095",
    #Ldh
        "23444",

    ]
    
    #Lqo 
    model.reactions["Lqo"] = Reaction(;
        name = "Lactate-quinone oxidoreductase",
        stoichiometry = Dict(
            "CHEBI:16374" => -1, #a menaquinone 
            "CHEBI:16651" => -1, #(S)-lactate
            "CHEBI:18151" => 1, #a menaquinol 
            "CHEBI:15361" => 1, #pyruvate
        ),
        annotations = Dict(
            "REACTION" => [
                "a menaquinone + (S)-lactate = a menaquinol + pyruvate"
            ],
        )
    )

    #cyt aa3 
    model.reactions["cyt_aa3"] = Reaction(;
        name = "Cytochrome aa3 oxidase",
        stoichiometry = Dict(
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -1, #H+
            "CHEBI:29033" => -1, #Fe(2+)
            "CHEBI:29034" => 1, #Fe(3+)
            "CHEBI:15377" => 1, #H2O
        ),
        annotations = Dict(
            "REACTION" => [
                "O2 + Fe(2+) + H+ = H2O + Fe(3+)"
            ]
        )
    )

    model.reactions["cyt_bd"] = Reaction(;
        name = "Cytochrome bd oxidase",
        stoichiometry = Dict(
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -1, #H+
            "CHEBI:18151" => -1, #a menaquinol
            "CHEBI:16374" => 1, #a menaquinone
            "CHEBI:15377" => 1, #H2O
            "CHEBI:15378_p" => 1, #H+ periplasm
        )
    )

    model.reactions["cyt_bo3"] = Reaction(;
        name = "Cytochrome bo3 oxidase",
        stoichiometry = Dict(
            "CHEBI:15379" => -1, #O2
            "CHEBI:15378" => -1, #H+
            "CHEBI:18151" => -1, #a menaquinol
            "CHEBI:16374" => 1, #a menaquinone
            "CHEBI:15377" => 1, #H2O
            "CHEBI:15378_p" => 1, #H+ periplasm
        )
    )



end
