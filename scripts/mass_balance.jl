using StaphylococcusAureus
import AbstractFBCModels as A
import ConstraintTrees as C
using COBREXA, HiGHS

model, reaction_isozymes = build_model();
model.reactions["EX_47013"].upper_bound = 0 #block ribose exchange
model.reactions["EX_15903"].upper_bound = 10 #limit glucose
sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)


model.reactions["EX_15379"].upper_bound = 100
sol = parsimonious_flux_balance_analysis(model;optimizer=HiGHS.Optimizer)

C.pretty(
    C.ifilter_leaves(sol.fluxes) do ix, x
        abs(x) > 1e-6 && startswith(string(last(ix)), "EX_")    
    end; 
    format_label = x -> A.reaction_name(model, string(last(x))),
)

Dict(x => y for (x,y) in sol.fluxes if abs(y)>100)



model.reactions["EX_15379"].upper_bound = 0


rids = filter(x -> !startswith(x, "EX_") && x != "biomass", A.reactions(model))
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
[rid for rid in unbal_rids if haskey(A.reaction_stoichiometry(model,rid),"CHEBI:29950")]

["20345","28037"]




for rid in unbal_rids 
    lb = copy(model.reactions[rid].lower_bound)
    ub = copy(model.reactions[rid].upper_bound)
    model.reactions[rid].lower_bound, model.reactions[rid].upper_bound = 0, 0
    sol = flux_balance_analysis(model; optimizer = HiGHS.Optimizer)
    if isnothing(sol) || sol.objective < 1e-3 
        model.reactions[rid].lower_bound, model.reactions[rid].upper_bound = lb, ub 
                println(rid)
    else 
        delete!(model.reactions,rid)
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

get_rid = ["15184","20213","16413","28767"]

model.reactions["28767"].stoichiometry["CHEBI:15379"] = 0.5

# macrochemical equation
eqn = Dict(
    "CHEBI:"*split(string(x),"_")[2] => -y for (x,y) in sol.fluxes if startswith(string(x),"EX") && abs(y)>1e-5
)
for (k,v) in A.reaction_stoichiometry(model,"biomass")
    if haskey(eqn,k)
        eqn[k] -= v * sol.objective 
    else
        eqn[k] = -v * sol.objective
    end
end
eqn
# mass balance
m = Dict()
for (k, v) in eqn
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
