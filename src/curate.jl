function curate!(model)

    # modify rhea reactions to use beta-D isomer instead of D-glucose 

    for (r,rxn) in model.reactions 
        if haskey(rxn.stoichiometry,"CHEBI:4167") # D-glucose -> beta-D-glucose
            rxn.stoichiometry["CHEBI:15903"] = rxn.stoichiometry["CHEBI:4167"]
            delete!(rxn.stoichiometry, "CHEBI:4167")
        end
        if haskey(rxn.stoichiometry,"CHEBI:61548") # D-glucose 6-phosphate -> β-D-glucose 6-phosphate 
            rxn.stoichiometry["CHEBI:58247"] = rxn.stoichiometry["CHEBI:61548"]
            delete!(rxn.stoichiometry,"CHEBI:61548")
        end
    end

    delete!(model.metabolites, "CHEBI:4167") # D-glucose
    delete!(model.metabolites, "CHEBI:61548") # D-glucose 6-phosphate

    # allow bidirectional H2O 
    model.reactions["EX_15377"].lower_bound = -1000
    model.reactions["EX_15377"].upper_bound = 1000
    model 
end
