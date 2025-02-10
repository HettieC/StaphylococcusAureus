using CSV, DataFrames, DataFramesMeta
using COBREXA, AbstractFBCModels
import AbstractFBCModels.CanonicalModel as CM
import COBREXA as X
using Literate
using RheaReactions

include("src/utils.jl")
model = CM.Model()

df = DataFrame(CSV.File("data/model/metabolic_reactions.csv"))


heteros = @rsubset(df, !ismissing(:Subunit))

gheteros = groupby(heteros, [:RHEA_ID, :Subunit])

gs = String[]
ms = RheaReactions.RheaMetabolite[]

dfs = gheteros #######

df = dfs[1]
#for df in dfs

    rid = parse(Int64,split(first(df.RHEA_ID),':')[2])
    grr = String.(df.Protein[:])
    stoich = Int.(df.Stoichiometry[:])
    append!(gs, grr)

    iso = X.Isozyme(; gene_product_stoichiometry=Dict(grr .=> stoich))

    if haskey(model.reactions, string(rid)) # isozyme

        push!(model.reactions[string(rid)].gene_association, iso)

    else # first time seeing this reaction

        rxn = get_reaction(rid)

        coeff_mets = get_reaction_metabolites(rid)
        stoichiometry = Dict(
            string(v.accession) => s
            for (s, v) in coeff_mets
        )

        append!(ms, last.(coeff_mets))

        ecs = isnothing(rxn.ec) ? df.EC : rxn.ec
        name = isnothing(rxn.name) ? df.reaction_name : rxn.name

        # direction
        # reversibility_index_threshold = 3.0
        # rev_ind = ismissing(first(df.RevIndex)) ? nothing : first(df.RevIndex)
        # dg = ismissing(first(df.DeltaG)) ? nothing : first(df.DeltaG)

        # if isnothing(rev_ind) || (abs(rev_ind) <= reversibility_index_threshold)
        #     lb = -1000
        #     ub = 1000
        # elseif rev_ind < -reversibility_index_threshold # forward
        #     lb = 0
        #     ub = 1000
        # elseif rev_ind > reversibility_index_threshold # reverse
        #     lb = -1000
        #     ub = 0
        # end

        model.reactions[string(rid)] = Reaction(;
            name=name,
            lower_bound=-1000.0,
            upper_bound=1000.0,
            #dg=dg,
            gene_association=[iso],
            stoichiometry=stoichiometry,
            annotations=Dict(
                "REACTION" => [rxn.equation],
                "EC" => ecs,
            ),
        )
    end
#end


add_genes!(model, gs)
add_metabolites!(model, ms)
