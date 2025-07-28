using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using HiGHS, JSON, CSV
using JSONFBCModels, DataFrames
using Latexify, CairoMakie

model, reaction_isozymes = build_model();
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)

carbon_sources = Dict(
    "15903" => "beta-D-glucose",
    "47013" => "D-ribose",
    "57972" => "alanine",
    "30031" => "succinate",
    "17754" => "glycerol",
    "29985" => "glutamate",
    "47013" => "ribose",
    "30089" => "acetate",
    "58723" => "glucosamine",
    "506227" => "n-acetyl-d-glucosamine",
    "18391" => "gluconate",
    "15589" => "malate",
    "29806" => "fumarate",
    "40886" => "L-arabinose",
    "15361" => "pyruvate"
)

model.reactions["EX_15903"].upper_bound = 0
model.reactions["EX_47013"].upper_bound = 0
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 200.0),
    ("membrane", membrane_gids, 130.0)
];

fba_growth = Dict{String,Float64}()
for (x,y) in carbon_sources 
    "CHEBI:$x" ∉ A.metabolites(model) && continue

    if !haskey(model.reactions,"EX_$x")
        println("rxn added: $y\n")
        model.reactions["EX_$x"] = CM.Reaction(;
            name = "$y exchange into cytosol",
            stoichiometry = Dict("CHEBI:$x" => 1),
            lower_bound = 0.0
        )
    end
    model.reactions["EX_$x"].upper_bound = 1/A.metabolite_formula(model,"CHEBI:$x")["C"]
    sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    fba_growth[y] = isnothing(sol) ? 0 : sol.objective 
    fba_growth[y] == 0 && continue

    model.reactions["EX_$x"].upper_bound = 0

end
fba_growth


ecfba_growth = Dict{String,Float64}()
acetate = Dict{String,Vector{Vector{Float64}}}()
for (x,y) in carbon_sources 
    "CHEBI:$x" ∉ A.metabolites(model) && continue

    model.reactions["EX_$x"].upper_bound = 1000
    ec_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    ecfba_growth[y] = isnothing(ec_sol) ? 0 : ec_sol.objective 
    isnothing(ec_sol) && continue 

    ### run overflow scan 
    model.reactions["biomass"].objective_coefficient = 0
    ac_flux = Float64[]
    vols = 0.1:copy(ec_sol.objective)/15:copy(ec_sol.objective)
    for biomass in vols
        model.reactions["biomass"].upper_bound = biomass
        model.reactions["biomass"].lower_bound = biomass-0.1

        ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity)
        
        ec_sol = optimized_values(
            ct;
            optimizer = HiGHS.Optimizer,
            objective = sum_value(ct.gene_product_amounts),
            sense = Minimal
        )
        push!(ac_flux,ec_sol.fluxes["EX_30089"])
    end
    acetate[y] = [vols,ac_flux]
    model.reactions["EX_$x"].upper_bound = 0
    model.reactions["biomass"].upper_bound = 1000 
    model.reactions["biomass"].lower_bound = 0
    model.reactions["biomass"].objective_coefficient = 1

end
ecfba_growth
acetate 

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

f = Figure(; size=(10cm, 6cm))

ax = Axis(
    f[1,1];
    backgroundcolor=:transparent,
    xlabel = "Growth rate (gDW/h)",
    ylabel = "Acetate exchange rate (mmol/gDW/h)",
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    #xticks = [0,0.5,1,1.5,2,2.5],
    ygridvisible=false,
    xgridvisible=false,
)
for (k,v) in acetate 
    k == "acetate" && continue
    lines!(
        ax,
        v[1],
        abs.(v[2]),
        label = k
    )
end

f

df = DataFrame(Carbon_source=String[],Growth=Float64[])
for (x,y) in fba_growth
    push!(
            df,
            [
                x,
                round(y,sigdigits=4),
            ])
end
df
sort!(df,:Growth)
latexify(df; env = :table, booktabs = true, latex = false) |> print

