function scale_biomass!(model)
    chebi_mass = Dict{String,Float64}()
    open("data/databases/chebi/chebi_core.obo","r") do io 
        i = 0
        chebi = "" 
        mass = 0
        for ln in eachline(io)
            i += 1
            i < 20 && continue
            if startswith(ln,"id: ")
                chebi = string(split(ln,"id: ")[2])
            elseif startswith(ln,"property_value: http://purl.obolibrary.org/obo/chebi/mass ")
                mass = parse(Float64,split(ln)[3][2:end-1])    
            elseif ln == "[Term]"
                chebi_mass[chebi] = mass
            end
        end
    end
    mass = 0
    for (m,s) in model.reactions["biomass"].stoichiometry 
        s > 0 && continue 
        mass += chebi_mass[m]*s
    end
    # mass in Da=g/mol, i want 1g/mmol so *1000
    model.reactions["biomass"].stoichiometry = Dict(x => -1000*y/mass for (x,y) in model.reactions["biomass"].stoichiometry)
    return model 
end
export scale_biomass!
