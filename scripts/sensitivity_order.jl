
flux_idxs = findall(x -> first(x) == :fluxes, vids)
flux_ids = last.(vids[flux_idxs])

# put flux ids into dictionary with their capacity allocation based on the gene_association_dnf of the reaction of the flux id
param_location = Dict{Symbol, String}()
for p_id in parameters
    if !isnothing(pruned_model.reactions[string(p_id)].gene_association_dnf)
        gids = vcat(pruned_model.reactions[string(p_id)].gene_association_dnf...)
        if gids == ["gid"]
            param_location[p_id] = "membrane"
        else
            param_location[p_id] = any(gid -> gid in membrane_gids, gids) ? "membrane" : "cytosol"
        end
    else
        param_location[p_id] = "unknown"
    end
end

# make heatmap of sens[flux_idxs,:]' with the flux_ids ordered by their flux_location and the param_ids in the same order as flux_ids

sort_perm = sortperm([param_location[p_id] for p_id in parameters])
# order the x-axis the same as the y-axis
sort_perm_dict = Dict(zip(parameters, sort_perm))

parameters[sort_perm]
# put the flux_ids in the same order as parameters[sort_perm]
order = Dict{Symbol, Int64}()
i = 0
for flux_id in flux_ids
    if flux_id ∈ parameters
        order[flux_id] = findfirst(x -> x == flux_id, parameters[sort_perm])
    else
        i += 1
        order[flux_id] = length(parameters) + i # put it at the end
    end
end
order

# sort the flux_ids by their order
sort_perm_fluxes = sortperm([order[flux_id] for flux_id in flux_ids])


# make yticklabels the flux_location of the flux_ids
yticklocs = [haskey(param_location,flux_id) ? param_location[flux_id] : "" for flux_id in flux_ids[sort_perm_fluxes]]
xticklocs = [param_location[p_id] for p_id in parameters[sort_perm]]


xticks = [
    0,
    findlast(x -> x == "cytosol", xticklocs)/2,
    findlast(x -> x == "cytosol", xticklocs),
    (length(parameters)+findlast(x -> x == "cytosol", xticklocs))/2,
    length(parameters),
]
xticklabels = ["0","cytosol","$(findlast(x -> x == "cytosol", xticklocs))","memb.",""]
yticks = [
    0,
    findlast(x -> x == "cytosol", yticklocs)/2,
    findlast(x -> x == "cytosol", yticklocs),
    (findlast(x -> x == "membrane", yticklocs)+findlast(x -> x == "cytosol", yticklocs))/2,
    findlast(x -> x == "membrane", yticklocs),
    (length(flux_ids)+findlast(x -> x == "membrane", yticklocs))/2,
    length(flux_ids),
]
yticklabels = ["0","cytosol","$(findlast(x -> x == "cytosol", yticklocs))","memb.","$(findlast(x -> x == "membrane", yticklocs))","NA",""]

#scale the sensitivities to be between -1 and 1
sens_scaled = sens[flux_idxs[sort_perm_fluxes],sort_perm]'  


f = Figure(; size=(12cm,10cm))#, backgroundcolor=:transparent)
ax = Axis(
    f[1,1],
    xlabel = L"\text{Enzyme, }p",
    yticklabelrotation = pi / 2,
    ylabel = L"Reaction flux sensitivity, $\frac{\partial v}{\partial p}$",   
    xticks = (xticks, xticklabels),
    xlabelsize=7pt,
    ylabelsize=7pt,
    xticklabelsize=6pt,
    yticklabelsize=6pt,
    yticks = (yticks, yticklabels),
)

#plot heatmap of sens[flux_idxs[sortperm(a)],:]' in the y-axis order of the sort_perm and

hm = heatmap!(
    ax,
    sens[flux_idxs[sort_perm_fluxes],sort_perm]',;
    colormap = reverse(ColorSchemes.RdBu),
    #colorrange = (-0.2,0.2)
)
# add bracket to the heatmap from start of x cytosol ticks to end of cytosol ticks

Colorbar(f[1, 2], hm, ticklabelsize = 5pt)
f

save("data/plots/whole_sens.png", f, px_per_unit = 1200/inch)
