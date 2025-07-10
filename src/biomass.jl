function scale_biomass!(model)

    total_mass = 0
    for (m,s) in model.reactions["biomass"].stoichiometry 
        total_mass -= parse(Float64, first(model.metabolites[m].annotations["molarmass"]))*s
    end

    # mass in Da=g/mol, i want 1g/mmol so *1000
    model.reactions["biomass"].stoichiometry = Dict(x => 1000*y/total_mass for (x,y) in model.reactions["biomass"].stoichiometry)
    return model 
end
export scale_biomass!
