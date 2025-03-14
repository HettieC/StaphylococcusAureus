
open("data/databases/ST398.txt") do io
    open("data/sequence.txt","w") do io2
        for ln in eachline(io)
            if startswith(ln, '>')
                println(io2,">$(split(split(ln, "locus_tag=")[2], ']'; limit=2)[1])")
            else
                println(io2,ln)
            end
        end
    end
end
