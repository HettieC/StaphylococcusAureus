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

growth
ac_flux 
bounds

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=20)
colors = Makie.wong_colors()[[2,1,6,4]]
colors = Dict("100" => colors[1],"110"=>colors[2],"120"=>colors[3],"130"=>colors[4])

xticks_enz = [0.0,1]
xticks_ac = [0.0,1]
xticks_glc = [0.0,5]
idxs = Int64[]
for (x,y) in ac_flux 
    if startswith(x,"enzyme")
        push!(xticks_enz,round(growth[x][findfirst(a -> a>1, y)],sigdigits=2))
    elseif startswith(x,"acetate")
        push!(xticks_ac,round(growth[x][findfirst(a -> a>0.1, y)],sigdigits=2))
    else 
        push!(xticks_glc,round(glucose[x][findfirst(a -> a>0.1, y)],sigdigits=2))
        push!(idxs,findfirst(a -> a>0.1, y))
    end
end
sort!(idxs)
xticklabels_glc = [
    "0",
    "5",
    "10",
    "$(round(glucose["glc:100"][idxs[1]],sigdigits=3))\n$(round(growth["glc:100"][idxs[1]],sigdigits=3))",
    "$(round(glucose["glc:110"][idxs[2]],sigdigits=3))\n$(round(growth["glc:110"][idxs[2]],sigdigits=3))",
    "$(round(glucose["glc:120"][idxs[3]],sigdigits=3))\n$(round(growth["glc:120"][idxs[3]],sigdigits=3)) gDW/h",
    "$(round(glucose["glc:130"][idxs[3]],sigdigits=3))\n$(round(growth["glc:130"][idxs[3]],sigdigits=3)) gDW/h"]

xticklabels_glc
fig_enz = Figure(; size=(9cm, 7.2cm));
ax_enz = Axis(
    fig_enz[1,1];
    backgroundcolor=:transparent,
    ylabel = "Acetate exchange rate (mmol/gDW/h)",
    xlabel = "Growth rate (gDW/h)",
    xlabelsize=10pt,
    ylabelsize=10pt,
    xticklabelsize=8pt,
    yticklabelsize=8pt,
    xticks=xticks_enz,
    ygridvisible=false,
    xgridvisible=false,
);
fig_ac = Figure(; size=(9cm, 7.2cm));
ax_ac = Axis(
    fig_ac[1,1];
    backgroundcolor=:transparent,
    ylabel = "Acetate exchange rate (mmol/gDW/h)",
    xlabel = "Growth rate (gDW/h)",
    xlabelsize=10pt,
    ylabelsize=10pt,
    xticklabelsize=8pt,
    yticklabelsize=8pt,
    xticks=xticks_ac,
    ygridvisible=false,
    xgridvisible=false,
);
fig_glc = Figure(; size=(9cm, 7.2cm));
ax_glc = Axis(
    fig_glc[1,1];
    backgroundcolor=:transparent,
    ylabel = "Acetate secretion (mmol/gDW/h)",
    xlabel = "\nGlucose uptake (mmol/gDW/h)",
    xlabelsize=10pt,
    ylabelsize=10pt,
    xticklabelsize=8pt,
    yticklabelsize=8pt,
    xticks = xticks_glc,
    #xticklabelrotation = pi/3,
    ygridvisible=false,
    xgridvisible=false,
);
for (x,y) in ac_flux 
    if startswith(x,"enzyme")
        lines!(
            ax_enz,
            growth[x],
            y,
            label = "$(split(x,":")[2])mg/gDW",
            linewidth=2.5,
            color = colors[split(x,":")[2]],
        )
        band!(
            ax_enz,
            growth[x][findfirst(((a,b),) -> a>0.999*parse(Float64,split(x,":")[2]),bounds[x])]:0.0001:growth[x][end],
            0,
            maximum(vcat(values(filter(((k,v),) -> startswith(k,"enz"),ac_flux))...)),
            color = (colors[split(x,":")[2]],0.2)
        )
    elseif startswith(x,"acetate") 
        lines!(
            ax_ac,
            growth[x],
            y,
            label = "$(split(x,":")[2])mg/gDW",
            linewidth=2.5,
            color = colors[split(x,":")[2]]
        )
        band!(
            ax_ac,
            growth[x][findfirst(((a,b),) -> a>0.999*parse(Float64,split(x,":")[2]),bounds[x])]:0.0001:growth[x][end],
            0,
            maximum(vcat(values(filter(((k,v),) -> startswith(k,"ace"),ac_flux))...)),
            color = (colors[split(x,":")[2]],0.2)
        )
    else 
        lines!(
            ax_glc,
            glucose[x],
            y,
            label = "$(split(x,":")[2])mg/gDW",
            linewidth=2.5,
            color = colors[split(x,":")[2]]
        )
        band!(
            ax_glc,
            glucose[x][findfirst(((a,b),) -> a>0.999*parse(Float64,split(x,":")[2]),bounds[x])]:0.0001:glucose[x][end],
            0,
            maximum(vcat(values(filter(((k,v),) -> startswith(k,"glc"),ac_flux))...)),
            color = (colors[split(x,":")[2]],0.2)
        )
    end
end
xlims!(ax_enz,(0,maximum(growth["enzymes:130"])))
ylims!(ax_enz, (0,maximum(vcat(values(filter(((k,v),) -> startswith(k,"enz"),ac_flux))...))))
xlims!(ax_ac,(0,maximum(growth["acetate:130"])))
ylims!(ax_ac, (0,maximum(vcat(values(filter(((k,v),) -> startswith(k,"ace"),ac_flux))...))))

fig_glc 

ax_glc_2 = Axis(
    fig_glc[1,1];
    backgroundcolor=:transparent,
    xlabelsize=6pt,
    xticklabelcolor=Makie.wong_colors()[3],
    ylabelsize=0pt,
    xticklabelsize=7pt,
    yticklabelsize=0pt,
    yticksvisible=false,
    xticksvisible=false,
    xticks=([7,glucose["glc:100"][idxs[1]],glucose["glc:110"][idxs[2]],glucose["glc:120"][idxs[3]],glucose["glc:130"][idxs[4]]],["\n\nGrowth (gDW/h):","\n\n$(round(growth["glc:100"][idxs[1]],sigdigits=2))","\n\n$(round(growth["glc:110"][idxs[2]],sigdigits=2))","\n\n$(round(growth["glc:120"][idxs[3]],sigdigits=2))","\n\n$(round(growth["glc:130"][idxs[4]],sigdigits=2))"]),
    ygridvisible=false,
    xgridvisible=false,
)
linkaxes!(ax_glc,ax_glc_2)
hidespines!(ax_glc_2)
xlims!(ax_glc_2,(0,maximum(glucose["glc:130"])))
ylims!(ax_glc_2, (0,maximum(vcat(values(filter(((k,v),) -> startswith(k,"glc"),ac_flux))...))))
labels = ["100 mg/gDW","110 mg/gDW","120 mg/gDW","130 mg/gDW"]
elements = [
    LineElement(color=colors["100"],linestyle=:solid,linewidth=3),
    LineElement(color=colors["110"],linestyle=:solid,linewidth=3),
    LineElement(color=colors["120"],linestyle=:solid,linewidth=3),
    LineElement(color=colors["130"],linestyle=:solid,linewidth=3)]
Legend(fig_glc[1,1],
elements,
labels,
labelsize=8pt,
tellheight = false,
tellwidth = false,
margin = (10, 10, 10, 10),
halign = :left, valign = :top,
)
fig_glc
fig_ac
fig_enz

save("data/plots/acetate_fix_glc.png",fig_glc,px_per_unit = 1200/inch)
save("data/plots/acetate_min_enz.png",fig_enz,px_per_unit = 1200/inch)
save("data/plots/acetate_min_ac.png",fig_ac,px_per_unit = 1200/inch)


set_theme!(figure_padding=5)
# plot acetate min and enzyme min on same plot 
fig = Figure(; size=(9cm, 7.2cm))
ax = Axis(
    fig[1,1];
    backgroundcolor=:transparent,
    ylabel = "Acetate exchange rate (mmol/gDW/h)",
    xlabel = "Growth rate (gDW/h)",
    xlabelsize=10pt,
    ylabelsize=10pt,
    xticklabelsize=8pt,
    yticklabelsize=8pt,
    xticks=xticks_enz,
    ygridvisible=false,
    xgridvisible=false,
);
for (x,y) in ac_flux 
    if startswith(x,"enzyme")
        lines!(
            ax,
            growth[x],
            y,
            label = "$x mg/gDW",
            linewidth=3,
            color = colors[split(x,":")[2]],
            linestyle = :dot,
        )
        band!(
            ax,
            growth[x][findfirst(((a,b),) -> a>0.999*parse(Float64,split(x,":")[2]),bounds[x])]:0.0001:growth[x][end],
            0,
            maximum(vcat(values(filter(((k,v),) -> startswith(k,"enz"),ac_flux))...)),
            color = (colors[split(x,":")[2]],0.2)
        )
    elseif startswith(x,"acetate") 
        lines!(
            ax,
            growth[x],
            y,
            label = "$x mg/gDW",
            linewidth=3,
            color = colors[split(x,":")[2]],
        )
    end
end
xlims!(ax,(0,maximum(growth["acetate:130"])*1.01))
ylims!(ax, (0,maximum(vcat(values(filter(((k,v),) -> !startswith(k,"glc"),ac_flux))...))))
labels = ["enzyme:100 mg/gDW","110 mg/gDW","120 mg/gDW","130 mg/gDW","acetate:100 mg/gDW","110 mg/gDW","120 mg/gDW","130 mg/gDW"]
elements = [
    LineElement(color=colors["100"],linestyle=:dot,linewidth=3), 
    LineElement(color=colors["110"],linestyle=:dot,linewidth=3),
    LineElement(color=colors["120"],linestyle=:dot,linewidth=3),
    LineElement(color=colors["130"],linestyle=:dot,linewidth=3),
    LineElement(color=colors["100"],linestyle=:solid,linewidth=3),
    LineElement(color=colors["110"],linestyle=:solid,linewidth=3),
    LineElement(color=colors["120"],linestyle=:solid,linewidth=3),
    LineElement(color=colors["130"],linestyle=:solid,linewidth=3)]
Legend(fig[1,1],
elements,
labels,
labelsize=8pt,
tellheight = false,
tellwidth = false,
margin = (10, 10, 10, 10),
halign = :left, valign = :top,
)

fig

save("data/plots/acetate_enzyme_and_acetate.png",fig,px_per_unit = 1200/inch)



