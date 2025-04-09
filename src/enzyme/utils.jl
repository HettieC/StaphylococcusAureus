function get_gene_product_molar_mass(gids)
    AA_mass = Dict(
        'A' => 89,
        'R' => 174,
        'N' => 132,
        'D' => 133,
        'B' => 133,
        'C' => 121,
        'Q' => 146,
        'E' => 147,
        'Z' => 147,
        'G' => 75,
        'H' => 155,
        'I' => 131,
        'L' => 131,
        'K' => 146,
        'M' => 149,
        'F' => 165,
        'P' => 115,
        'S' => 105,
        'T' => 119,
        'W' => 204,
        'Y' => 181,
        'V' => 117,
    )
    seqs = Dict{String,String}()
    open("data/sequence.txt","r") do io
        prev = ""
        for ln in eachline(io)
            if startswith(ln,'>')
                seqs[ln[2:end]] = ""
                prev = ln[2:end]
            else 
                seqs[prev] = seqs[prev]*ln 
            end
        end
    end

    gene_product_molar_mass = Dict{String,Float64}()
    for gid in gids 
        gene_product_molar_mass[gid] = sum([AA_mass[aa] for aa in seqs[gid]])
    end

    println(110*length(seqs["SAPIG1968"]))
    return gene_product_molar_mass     
    
end

