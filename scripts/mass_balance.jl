using StaphylococcusAureus
import AbstractFBCModels as A
import ConstraintTrees as C
using COBREXA

model, reaction_isozymes = build_model()

ids = filter(x -> !startswith(x, "EX_") && x != "biomass", A.reactions(model))
unbal_rids = String[]
for rid in rids
    s = A.reaction_stoichiometry(model, rid)
    m = Dict()
    for (k, v) in s
        if isnothing(A.metabolite_formula(model,k))
            println(rid)
            println(k)
            continue 
        end
        for (kk, vv) in A.metabolite_formula(model, k)
            m[kk] = get(m, kk, 0) + vv * v
        end
    end
    m
    all(values(m) .== 0) || push!(unbal_rids, rid) 
end
unbal_rids
[rid for (rid,v) in sol.fluxes if abs(v)>1e-5 && string(rid) ∈ unbal_rids]

for rid in unbal_rids 
    lb = copy(model.reactions[rid].lower_bound)
    ub = copy(model.reactions[rid].upper_bound)
    model.reactions[rid].lower_bound, model.reactions[rid].upper_bound = 0, 0
    sol = flux_balance_analysis(model; optimizer = HiGHS.Optimizer)
    if isnothing(sol) || sol.objective < 1e-3 
        model.reactions[rid].lower_bound, model.reactions[rid].upper_bound = lb, ub 
    else 
        delete!(model.reactions,rid)
        println(rid)
    end    
end

sol = flux_balance_analysis(model; optimizer = HiGHS.Optimizer)

unbal_stoich = Dict{String,Dict{String,Float64}}()
for rid in unbal_rids
    s = A.reaction_stoichiometry(model, rid)
    m = Dict{String,Float64}()
    for (k, v) in s
        if isnothing(A.metabolite_formula(model,k))
            println(rid)
            println(k)
            continue 
        end
        for (kk, vv) in A.metabolite_formula(model, k)
            m[kk] = get(m, kk, 0) + vv * v
        end
    end
    m
    unbal_stoich[rid] = m
end
unbal_stoich



rid = "cyt_bd"
s = A.reaction_stoichiometry(model, rid)
m = Dict()
for (k, v) in s
    if isnothing(A.metabolite_formula(model,k))
        println(rid)
        println(k)
        continue 
    end
    for (kk, vv) in A.metabolite_formula(model, k)
        m[kk] = get(m, kk, 0) + vv * v
    end
end
m
