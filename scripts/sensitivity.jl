using StaphylococcusAureus
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import ConstraintTrees as C
using COBREXA
import COBREXA as X
using CairoMakie
using HiGHS, JSON, CSV
using JSONFBCModels, DataFrames
import DifferentiableMetabolism as D
import FastDifferentiation as F
const Ex = F.Node
using CairoMakie
using Latexify
flux_zero_tol = 1e-6
gene_zero_tol = 1e-6

model, reaction_isozymes = build_model()

gene_product_molar_masses, membrane_gids = enzyme_constraints!(model,reaction_isozymes)

escher_model = change_reaction_names(model)
save_model(convert(JSONFBCModels.JSONFBCModel, escher_model), "data/escher_model.json")
model.reactions["EX_16236"].lower_bound = 0 #block ethanol exchange
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_16651"].lower_bound = 0 #block (s)-lactate exchange
model.reactions["EX_16004"].lower_bound = 0 #block (r)-lactate exchange
model.reactions["EX_15740"].lower_bound = 0 #block formate exchange
model.reactions["EX_15378"].lower_bound = 0 #block H+ exchange


capacity = [
    ("cytosol", [g for g in A.genes(model) if g ∉ membrane_gids], 300.0),
    ("membrane", membrane_gids, 120.0)
];

ec_sol = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)

# atp-producing rxns
C.pretty(
    C.ifilter_leaves(ec_sol.fluxes) do ix, x
        abs(x) > 1e-5 && 
            haskey(model.reactions[string(last(ix))].stoichiometry,"CHEBI:30616") && 
            ((model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] > 0 && x > 1e-5) || 
            (model.reactions[string(last(ix))].stoichiometry["CHEBI:30616"] < 0 && x < -1e-5))
    end; 
    format_label = x -> (string(last(x)),A.reaction_name(model, string(last(x)))),
)



pruned_model, pruned_reaction_isozymes = D.prune_model(
    model,
    ec_sol.fluxes,
    ec_sol.gene_product_amounts,
    reaction_isozymes,
    ec_sol.isozyme_forward_amounts,
    ec_sol.isozyme_reverse_amounts,
    flux_zero_tol,
    gene_zero_tol,
);

pruned_sol = enzyme_constrained_flux_balance_analysis(
    pruned_model;
    reaction_isozymes = pruned_reaction_isozymes,
    gene_product_molar_masses,
    capacity,
    optimizer=HiGHS.Optimizer,
)

Dict(x=>[y,abs(ec_sol.fluxes[x])] for (x,y) in pruned_sol.fluxes if abs(abs(y)-abs(ec_sol.fluxes[x]))>1e-5)


# differentiation variables
parameter_values = Dict(
    Symbol(x) => y[iso_id].kcat_forward for (x,y) in pruned_reaction_isozymes for (iso_id,iso) in y if ec_sol.isozyme_forward_amounts[x][iso_id] > 1e-8 || ec_sol.isozyme_reverse_amounts[x][iso_id] > 1e-8
)

parameters = Symbol.(collect(keys(parameter_values)))
rid_pid = Dict(rid => Ex(Symbol(rid)) for (rid, v) in pruned_reaction_isozymes)
rid_gcounts = Dict(rid => [v.gene_product_stoichiometry for (k, v) in d][1] for (rid, d) in pruned_reaction_isozymes)

parameter_isozymes = Dict{String,Dict{String,X.IsozymeT{Ex}}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
parameter_isozymes = Dict(
    rid => Dict(
        iso_id => X.IsozymeT{Ex}(
            gene_product_stoichiometry = iso.gene_product_stoichiometry,
            kcat_forward = rid_pid[rid],
            kcat_reverse = nothing
        )
        for (iso_id,iso) in y
    ) 
    for (rid,y) in pruned_reaction_isozymes
)

pkm = X.enzyme_constrained_flux_balance_constraints( # pruned kinetic model
    pruned_model;
    reaction_isozymes = pruned_reaction_isozymes,
    gene_product_molar_masses,
    capacity
)

pruned_solution = D.optimized_values(
    pkm,
    parameter_values;
    objective = pkm.objective.value,
    optimizer = HiGHS.Optimizer,
)

pkm_kkt, vids = D.differentiate_prepare_kkt(pkm, pkm.objective.value, parameters)
sens = D.differentiate_solution(
    pkm_kkt,
    pruned_solution.primal_values,
    pruned_solution.equality_dual_values,
    pruned_solution.inequality_dual_values,
    parameter_values,
    scale = true, # unitless sensitivities
)


### use finite diff...
## loop over parameters and change kcats one by one
pruned_reaction_isozymes = Dict(r => Dict("iso" => iso for (id,iso) in isos) for (r,isos) in pruned_reaction_isozymes)
delta_pos = 1.001
delta_neg = 0.999

d_v = zeros(length(pruned_reaction_isozymes),length(A.reactions(pruned_model)))
for (i,r) in enumerate(collect(keys(pruned_reaction_isozymes)))
    println(r)
    pruned_reaction_isozymes[r]["iso"].kcat_forward *= delta_pos
    # solve the gecko model with new isozyme
    ec_solution_new = X.enzyme_constrained_flux_balance_analysis(
        pruned_model;
        reaction_isozymes=pruned_reaction_isozymes,
        gene_product_molar_masses=gene_product_molar_masses,
        capacity=capacity,
        optimizer = HiGHS.Optimizer
    )

    v_pos = collect(values(ec_solution_new.fluxes))

    pruned_reaction_isozymes[r]["iso"].kcat_forward *= delta_neg / delta_pos

    # solve the gecko model with new isozyme
    ec_solution_new = X.enzyme_constrained_flux_balance_analysis(
        pruned_model;
        reaction_isozymes=pruned_reaction_isozymes,
        gene_product_molar_masses=gene_product_molar_masses,
        capacity=capacity,
        optimizer=HiGHS.Optimizer
    )

    v_neg = collect(values(ec_solution_new.fluxes))


    #return to original value
    pruned_reaction_isozymes[r]["iso"].kcat_forward /= delta_neg

    # scaled sensitivity
    d_v[i,:] = pruned_reaction_isozymes[r]["iso"].kcat_forward*(
        v_pos .- v_neg
    ) ./ ((delta_pos - delta_neg) * pruned_sol.fluxes[r])
end

heatmap(d_v)

# pathways of reactions 
minimum(d_v)

df = DataFrame(d_v, :auto)
collect(keys(pruned_reaction_isozymes))
collect(values(ec_solution_new.fluxes))

rename!(df,Symbol.(collect(keys(pruned_sol.fluxes))))
oxphos_reactions = ["Sdh", "Mqo", "cyt_bo3",]
df[:,Symbol.(oxphos_reactions)]
heatmap(Matrix(df[:,Symbol.(oxphos_reactions)]))
