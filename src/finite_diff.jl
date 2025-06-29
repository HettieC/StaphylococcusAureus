"""
Use finite difference on a model.
"""
function finite_difference(pruned_model, pruned_reaction_isozymes, gene_product_molar_masses, capacity;flux_zero_tol = 1e-6,gene_zero_tol = 1e-6)
    using HiGHS

    ## loop over parameters and change kcats one by one
    delta_pos = 1.001
    delta_neg = 0.999

    d_v = zeros(length(pruned_reaction_isozymes),length(A.reactions(pruned_model)))
        for (i,r) in enumerate(collect(keys(pruned_reaction_isozymes)))
            i > 3 && break
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
            d_v[i,:] = (pruned_reaction_isozymes[r]["iso"]/pruned_sol.fluxes[r])*(
                v_pos .- v_neg
            ) ./ ((delta_pos - delta_neg) * pruned_reaction_isozymes[r]["iso"].kcat_forward)
        end
    return d_v
end
export finite_difference
