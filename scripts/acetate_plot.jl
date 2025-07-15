using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
using CairoMakie
using HiGHS, JSON
using JSONFBCModels

model, reaction_isozymes = build_model();
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)

model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange

# loop through 
ac_flux = Dict{String,Vector{Float64}}()
growth = Dict{String,Vector{Float64}}()
bounds = Dict{String,Vector{Tuple{Float64,Float64}}}()
for memb in [100,110,120]
    model.reactions["biomass"].objective_coefficient = 1
    model.reactions["biomass"].upper_bound = 1000
    capacity = [
        ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 200.0),
        ("membrane", membrane_gids, memb)
    ];
    first_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    model.reactions["biomass"].objective_coefficient = 0
    for obj in ["acetate","enzymes"]
        ac_flux["$obj:$memb"] = Float64[]
        growth["$obj:$memb"] = Float64[]
        bounds["$obj:$memb"] = Vector{Tuple{Float64,Float64}}()
        if obj == "enzymes"
            for bio in 0:0.01:copy(first_sol.objective)
                model.reactions["biomass"].upper_bound = bio
                model.reactions["biomass"].lower_bound = bio-0.02
                ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity)
                new_sol = optimized_values(
                    ct;
                    optimizer = HiGHS.Optimizer,
                    objective = sum_value(ct.gene_product_capacity),
                    sense = Minimal
                )
                push!(ac_flux["$obj:$memb"],abs(new_sol.fluxes["EX_30089"]))
                push!(bounds["$obj:$memb"],(new_sol.gene_product_capacity.membrane,new_sol.gene_product_capacity.cytosol))
                push!(growth["$obj:$memb"],new_sol.fluxes["biomass"])
            end
        else
            for bio in 0:0.01:copy(first_sol.objective)
                model.reactions["biomass"].upper_bound = bio
                model.reactions["biomass"].lower_bound = bio-0.02            
                ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity)
                new_sol = optimized_values(
                    ct;
                    optimizer = HiGHS.Optimizer,
                    objective = ct.fluxes.EX_30089.value,
                    sense = Maximal
                )
                push!(ac_flux["$obj:$memb"],abs(new_sol.fluxes["EX_30089"]))
                push!(bounds["$obj:$memb"],(new_sol.gene_product_capacity.membrane,new_sol.gene_product_capacity.cytosol))
                push!(growth["$obj:$memb"],new_sol.fluxes["biomass"])  
            end
        end
    end
end

ac_flux
glucose = Float64[]
for memb in [100,110,120]
    capacity = [
        ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 200.0),
        ("membrane", membrane_gids, memb)
    ];
    first_sol = enzyme_constrained_flux_balance_analysis(
        model;
        reaction_isozymes,
        gene_product_molar_masses,
        capacity,
        optimizer=HiGHS.Optimizer,
    )
    ac_flux["glc:$memb"] = Float64[]
    growth["glc:$memb"] = Float64[]
    bounds["glc:$memb"] = Vector{Tuple{Float64,Float64}}()
    for glc in 0:0.01:copy(first_sol.fluxes["EX_15903"])
        model.reactions["EX_15903"].upper_bound = glc
        model.reactions["EX_15903"].lower_bound = glc-0.02
        ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity)
        new_sol = optimized_values(
            ct;
            optimizer = HiGHS.Optimizer,
            objective = sum_value(ct.gene_product_capacity),
            sense = Minimal
        )
        isnothing(new_sol) && continue
        push!(ac_flux["glc:$memb"],abs(new_sol.fluxes["EX_30089"]))
        push!(bounds["glc:$memb"],(new_sol.gene_product_capacity.membrane,new_sol.gene_product_capacity.cytosol))
        push!(growth["glc:$memb"],new_sol.fluxes["biomass"])
        push!(glucose,new_sol.fluxes["EX_15903"])
    end
end

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)
colors = Makie.wong_colors()[[4,1,6]]

xtickstop = [0.0,1,2]
xticksbottom = [0.0,1,2]
for (x,y) in ac_flux 
    if startswith(x,"enzyme")
        push!(xtickstop,growth[x][findfirst(a -> a>1, y)])
    else
        push!(xticksbottom,growth[x][findfirst(a -> a>0.1, y)])
    end
end

f = Figure(; size=(12cm, 12cm))
ga = f[1, 1] = GridLayout()
axtop = Axis(
    ga[1,1];
    backgroundcolor=:transparent,
    ylabel = "Acetate exchange rate (mmol/gDW/h)",
    title = "Minimise enzyme usage",
    titlesize = 8pt,
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    xticks=xtickstop,
    ygridvisible=false,
    xgridvisible=false,
);
axbottom = Axis(
    ga[2,1];
    backgroundcolor=:transparent,
    title = "Minimise acetate secretion",
    xlabel = "Growth rate (gDW/h)",
    ylabel = "Acetate exchange rate (mmol/gDW/h)",
    titlesize=8pt,
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    xticks = xticksbottom,
    ygridvisible=false,
    xgridvisible=false,
);
ax3 = Axis(
    ga[1,3];
    backgroundcolor=:transparent,
    title = "Fix glucose, maximise biomass",
    xlabel = "Glucose uptake (mmol/gDW/h)",
    ylabel = "Acetate exchange rate (mmol/gDW/h)",
    titlesize=8pt,
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    xticks = xticksbottom,
    ygridvisible=false,
    xgridvisible=false,
)
xlims!(axbottom,(0,maximum(vcat(values(growth)...))))
ylims!(axbottom, (0,maximum(vcat(values(ac_flux)...))))
linkaxes!(axtop,axbottom)
for (x,y) in ac_flux 
    if startswith(x,"enzyme")
        lines!(
            axtop,
            growth[x],
            y,
            label = "Membrane: $(split(x,":")[2])mg/gDW",
            linewidth=2.5,
            color = split(x,":")[2] == "100" ? colors[1] : (split(x,":")[2] == "110" ? colors[2] : colors[3]),
        )
        band!(
            axtop,
            growth[x][findfirst(((a,b),) -> a>0.999*parse(Float64,split(x,":")[2]),bounds[x])]:0.0001:growth[x][end],
            0,
            maximum(vcat(values(ac_flux)...)),
            color = (split(x,":")[2] == "100" ? colors[1] : (split(x,":")[2] == "110" ? colors[2] : colors[3]),0.3)
        )
    else 
        lines!(
            axbottom,
            growth[x],
            y,
            label = "Membrane: $(split(x,":")[2])mg/gDW",
            linewidth=2.5,
            color = split(x,":")[2] == "100" ? colors[1] : (split(x,":")[2] == "110" ? colors[2] : colors[3])
        )
        band!(
            axbottom,
            growth[x][findfirst(((a,b),) -> a>0.999*parse(Float64,split(x,":")[2]),bounds[x])]:0.0001:growth[x][end],
            0,
            maximum(vcat(values(ac_flux)...)),
            color = (split(x,":")[2] == "100" ? colors[1] : (split(x,":")[2] == "110" ? colors[2] : colors[3]),0.3)
        )
    end
end

axislegend(
    axtop,
    position=:lt,
    labelsize = 5pt,
)
display(f)

save("data/plots/acetate_min_enzyme.png",f,px_per_unit = 1200/inch)

