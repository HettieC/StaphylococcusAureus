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

# loop through to colect data
ac_flux = Dict{String,Vector{Float64}}()
growth = Dict{String,Vector{Float64}}()
bounds = Dict{String,Vector{Tuple{Float64,Float64}}}()
for memb in [100,110,120,130]
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
            for bio in 0:0.1:copy(first_sol.objective)
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
            for bio in 0:0.1:copy(first_sol.objective)
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

model, reaction_isozymes = build_model();
gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange

glucose = Dict{String,Vector{Float64}}()
for memb in [100,110,120,130]
    model.reactions["EX_15903"].upper_bound = 1000 #unlimit glucose
    model.reactions["EX_15903"].lower_bound = 0 #unlimit glucose
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
    println(first_sol.fluxes["EX_15903"])
    ac_flux["glc:$memb"] = Float64[]
    growth["glc:$memb"] = Float64[]
    bounds["glc:$memb"] = Vector{Tuple{Float64,Float64}}()
    glucose["glc:$memb"] = Float64[]
    for glc in 0:0.05:copy(first_sol.fluxes["EX_15903"])
        model.reactions["EX_15903"].upper_bound = glc
        model.reactions["EX_15903"].lower_bound = glc-0.02
        ct = enzyme_constrained_flux_balance_constraints(model;reaction_isozymes,gene_product_molar_masses,capacity)
        new_sol = optimized_values(
            ct;
            optimizer = HiGHS.Optimizer,
            objective = ct.objective.value,
            sense = Maximal
        )
        isnothing(new_sol) && continue
        push!(ac_flux["glc:$memb"],abs(new_sol.fluxes["EX_30089"]))
        push!(bounds["glc:$memb"],(new_sol.gene_product_capacity.membrane,new_sol.gene_product_capacity.cytosol))
        push!(growth["glc:$memb"],new_sol.fluxes["biomass"])
        push!(glucose["glc:$memb"],new_sol.fluxes["EX_15903"])
    end
end

# make the plot 
inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=5)
colors = Makie.wong_colors()[[2,1,6,4]]
colors = Dict("100" => colors[1],"110"=>colors[2],"120"=>colors[3],"130"=>colors[4])

# figure of all three objectives
fig = Figure(; size=(12cm, 8cm));
ax = Axis(
    fig[1,1];
    backgroundcolor=:transparent,
    ylabel = "Acetate exchange rate (mmol/gDW/h)",
    xlabel = "Growth rate (gDW/h)",
    xlabelsize=10pt,
    ylabelsize=10pt,
    xticklabelsize=8pt,
    yticklabelsize=8pt,
    ygridvisible=false,
    xgridvisible=false,
);
for (x,y) in ac_flux 
    if startswith(x,"enzyme")
        
        lines!(
            ax,
            growth[x]/10,
            y,
            label = "$x mg/gDW",
            linewidth=3,
            color = colors[split(x,":")[2]],
            linestyle=:dash
        )
    elseif startswith(x,"acetate") 
        lines!(
            ax,
            growth[x]/10,
            y,
            label = "$x mg/gDW",
            linewidth=5,
            color = (colors[split(x,":")[2]],0.6),
            linestyle = :dot
        )

    elseif startswith(x,"glc")
        lines!(
            ax,
            growth[x]/10,
            y,
            label = "$x mg/gDW",
            linewidth=3,
            color = (colors[split(x,":")[2]],0.6),
            linestyle = :solid
        )
        band!(
            ax,
            growth[x][findfirst(((a,b),) -> a>0.999*parse(Float64,split(x,":")[2]),bounds[x])]/10:0.0001:growth[x][end]/10,
            0,
            maximum(vcat(values(filter(((k,v),) -> startswith(k,"enz"),ac_flux))...)),
            color = (colors[split(x,":")[2]],0.2),
        )
    end
end
fig
xlims!(ax,(0,maximum(growth["glc:130"]/10)*1.01))
ylims!(ax, (-0.1,maximum(vcat(values(filter(((k,v),) -> !startswith(k,"glc"),ac_flux))...))))
fig
labels = ["enzyme min","acetate min","glucose fix"]
elements = [
    LineElement(color=:gray,linestyle=(:dash,:dense),linewidth=3), 
    LineElement(color=:gray,linestyle=:dot,linewidth=3),
    LineElement(color=:gray,linestyle=:solid,linewidth=3),
]
Legend(fig[1,1],
    elements,
    labels,
    labelsize=8pt,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    halign = :left, valign = :top,
    "LP objective",
    titlesize=9
)
labels = ["100mg/DW","110mg/DW","120mg/gDW","130mg/gDW"]
elements = [
    PolyElement(color=colors["100"]), 
    PolyElement(color=colors["110"]),
    PolyElement(color=colors["120"]),
    PolyElement(color=colors["130"]), 
]
Legend(fig[1,1],
    elements,
    labels,
    labelsize=8pt,
    tellheight = false,
    tellwidth = false,
    margin = (10, 10, 10, 10),
    halign = :left, valign = :bottom,
    "Membrane bound",
    titlesize=9
)
fig

save("data/plots/overflow_scan.png",fig,px_per_unit = 1200/inch)


# separate figures 
objs = Dict("enzymes"=>:dash,"acetate"=>:dot,"glc"=>:solid) 
for (obj,ln) in objs
    fig = Figure(; size=(12cm, 8cm));
    ax = Axis(
        fig[1,1];
        backgroundcolor=:transparent,
        ylabel = "Acetate exchange rate (mmol/gDW/h)",
        xlabel = "Growth rate (gDW/h)",
        xlabelsize=10pt,
        ylabelsize=10pt,
        xticklabelsize=8pt,
        yticklabelsize=8pt,
        ygridvisible=false,
        xgridvisible=false,
    );
    for (x,y) in ac_flux 
        if startswith(x,obj)
            lines!(
                ax,
                growth[x]/10,
                y,
                label = "$x mg/gDW",
                linewidth=3,
                color = colors[split(x,":")[2]],
                linestyle=ln
            )
            band!(
                ax,
                growth[x][findfirst(((a,b),) -> a>0.999*parse(Float64,split(x,":")[2]),bounds[x])]/10:0.0001:growth[x][end]/10,
                0,
                maximum(vcat(values(filter(((k,v),) -> startswith(k,"enz"),ac_flux))...)),
                color = (colors[split(x,":")[2]],0.2),
            )
        end
    end
    xlims!(ax,(0,maximum(growth["$obj:130"]/10)*1.01))
    ylims!(ax, (-0.1,maximum(vcat(values(filter(((k,v),) -> !startswith(k,obj),ac_flux))...))))
    labels = ["100mg/DW","110mg/DW","120mg/gDW","130mg/gDW"]
    elements = [
        PolyElement(color=colors["100"]), 
        PolyElement(color=colors["110"]),
        PolyElement(color=colors["120"]),
        PolyElement(color=colors["130"]), 
    ]
    Legend(fig[1,1],
        elements,
        labels,
        labelsize=8pt,
        tellheight = false,
        tellwidth = false,
        margin = (10, 10, 10, 10),
        halign = :left, 
        valign = :center,
        "Membrane bound",
        titlesize=9
    )
    fig
    save("data/plots/overflow_scan_$obj.png",fig,px_per_unit = 1200/inch)
end
