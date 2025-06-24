"""
$(TYPEDSIGNATURES)
Get the reaction name, stoichiometry, and database cross references for a KEGG reaction ID
"""
function get_kegg_info(rxn_id::String;cache=true, max_retries = 5)
    rxn = nothing
    if contains(rxn_id, "(G)")#
        return nothing
    else
        if _is_cached("reaction",rxn_id)
            return _get_cache("reaction",rxn_id)
        else 
            retry_counter = 0 
            req = nothing
            while retry_counter <= max_retries
                retry_counter += 1
                try
                    req = HTTP.request("GET", "https://rest.kegg.jp/get/$rxn_id")
                catch
                    print("No KEGG entry matching this id: $rxn_id \n")
                    return nothing
                end
                out = Dict{String,Any}()
                lines = split(String(req.body), "\n")
                if split(lines[1])[3] != "Reaction"
                    throw(error("Entry $rxn_id not a reaction"))
                else
                    for ln in lines
                        if startswith(ln, "NAME")
                            out["name"] = strip(String(split(ln; limit = 2)[2]))
                        elseif startswith(ln, "PATHWAY")
                            out["pathway"] = [String(strip(split(ln; limit = 2)[2]))]
                        elseif startswith(strip(ln), "rn")
                            push!(out["pathway"], String(strip(ln)))
                        end
                    end
                end
                rxn = Dict(
                        "id" => rxn_id,
                        "name" => haskey(out, "name") ? out["name"] : nothing,
                        "pathway" => haskey(out, "pathway") ? out["pathway"] : nothing,
                    )
                if cache
                    _cache("reaction",rxn_id,rxn)
                end
            end
            return rxn
        end
    end
end
export get_kegg_info


"""
$(TYPEDSIGNATURES)
Clear the entire cache.
"""
clear_cache!() = begin
    for dir in readdir(CACHE_LOCATION)
        rm(joinpath(CACHE_LOCATION, dir), recursive = true)
        dir != "version.txt" && mkdir(joinpath(CACHE_LOCATION, dir)) # add back the empty dir
    end
    write(joinpath(CACHE_LOCATION, "version.txt"), string(Base.VERSION))    
    Term.tprint("{blue} Cache cleared! {/blue}")
end

"""
$(TYPEDSIGNATURES)
Checks if the reaction has been cached.
"""
_is_cached(database::String, id) = 
    isfile(joinpath(CACHE_LOCATION, database, string(id)))

"""
$(TYPEDSIGNATURES)
Return the cached reaction object.
"""
_get_cache(database::String, id) =
    deserialize(joinpath(CACHE_LOCATION, database, string(id)))

"""
$(TYPEDSIGNATURES)
Cache reaction object.
"""
_cache(database::String, id, item) =
    serialize(joinpath(CACHE_LOCATION, database, string(id)), item)
