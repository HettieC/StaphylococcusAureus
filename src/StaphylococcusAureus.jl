module StaphylococcusAureus

using DocStringExtensions
using CSV, DataFrames, DataFramesMeta
using COBREXA, AbstractFBCModels
import AbstractFBCModels as A
import AbstractFBCModels.CanonicalModel as CM
import COBREXA as X
using RheaReactions
using JSONFBCModels
using XLSX
using JSON
using HTTP, DocStringExtensions, Scratch, Serialization # all for kegg

# set up a cache for kegg reaction
const CACHE_DIRS = ["reaction", "reaction_metabolites"]
CACHE_LOCATION = ""

function __init__()
    global CACHE_LOCATION = @get_scratch!("kegg_data")

    for dir in CACHE_DIRS
        !isdir(joinpath(CACHE_LOCATION, dir)) && mkdir(joinpath(CACHE_LOCATION, dir))
    end

end

include("utils.jl")
include("annotate_utils.jl")
include("curate.jl")
include("reconstruct.jl")
include("specials.jl")
include("transporters.jl")
include("enzyme_constrain.jl")
include("kegg.jl")
include("finite_diff.jl")
include("biomass.jl")

end
