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
model.reactions["EX_16236"].lower_bound = 0 #block ethanol exchange
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_16651"].lower_bound = 0 #block (s)-lactate exchange
model.reactions["EX_16004"].lower_bound = 0 #block (r)-lactate exchange
model.reactions["EX_15740"].lower_bound = 0 #block formate exchange
model.reactions["EX_15378"].lower_bound = 0 #block H+ exchange

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 400.0),
    ("membrane", membrane_gids, 120.0)
];

fba_growth = Dict{String,Float64}()
ecfba_growth = Dict{String,Float64}()
acetate = Dict{String,Vector{Float64}}()
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
    model.reactions["EX_$x"].upper_bound = 10
    sol = flux_balance_analysis(model;optimizer=HiGHS.Optimizer)
    fba_growth[y] = isnothing(sol) ? 0 : sol.objective 
    delete!(model.reactions,"EX_$x")
    fba_growth[y] == 0 && continue

    # model.reactions["EX_$x"].upper_bound = 1000

    # ec_sol = enzyme_constrained_flux_balance_analysis(
    #     model;
    #     reaction_isozymes,
    #     gene_product_molar_masses,
    #     capacity,
    #     optimizer=HiGHS.Optimizer,
    # )
    # ecfba_growth[y] = isnothing(ec_sol) ? 0 : ec_sol.objective 
    # if abs(ecfba_growth[y]) < 1e-5 
    #     capacity = [
    #         ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 400.0),
    #         ("membrane", membrane_gids, 200.0)
    #     ];
    #     ec_sol = enzyme_constrained_flux_balance_analysis(
    #         model;
    #         reaction_isozymes,
    #         gene_product_molar_masses,
    #         capacity,
    #         optimizer=HiGHS.Optimizer,
    #     )
    #     ecfba_growth[y] = isnothing(ec_sol) ? 0 : ec_sol.objective 
    #     abs(ecfba_growth[y]) < 1e-5 && continue
    #     println("membrane 200: $y")
    # end

    ### run overflow scan 
    # acetate[y] = Float64[]
    # vols = 0.1:0.1:copy(ec_sol.objective)
    # for biomass in vols
    #     model.reactions["biomass"].upper_bound = biomass
    #     model.reactions["biomass"].lower_bound = biomass-0.1

    #     ec_sol = enzyme_constrained_flux_balance_analysis(
    #         model;
    #         reaction_isozymes,
    #         gene_product_molar_masses,
    #         capacity,
    #         optimizer=HiGHS.Optimizer,
    #     )
    #     push!(acetate[y],ec_sol.fluxes["EX_30089"])
    # end

    # f = Figure(; size=(10cm, 6cm))

    # ax = Axis(
    #     f[1,1];
    #     backgroundcolor=:transparent,
    #     xlabel = "Growth rate (gDW/h)",
    #     ylabel = "Acetate exchange rate (mMol/h)",
    #     title = y,
    #     xlabelsize=6pt,
    #     ylabelsize=6pt,
    #     xticklabelsize=5pt,
    #     yticklabelsize=5pt,
    #     ygridvisible=false,
    #     xgridvisible=false,
    # )
    # lines!(
    #     ax,
    #     vols,
    #     abs.(round.(acetate[y])),
    #     label = "Acetate exchange"
    # )
    # display(f)
    #model.reactions["biomass"].upper_bound = 1000
    # capacity = [
    #     ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 400.0),
    #     ("membrane", membrane_gids, 120.0)
    # ];
end
fba_growth
ecfba_growth

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

