using CSV, DataFrames, DataFramesMeta

#### make isozyme df 
iso = DataFrame(CSV.File("data/isozymes.csv"))
rxns = DataFrame(CSV.File("data/model/metabolic_reactions.csv"))
select!(rxns, :KEGG_ID, :RHEA_ID, :Protein)

df = leftjoin(rxns,iso, on = [:Protein, :KEGG_ID])
CSV.write("df.csv", df)

