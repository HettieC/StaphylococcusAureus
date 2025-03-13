using CSV, DataFrames, CairoMakie

df = DataFrame(CSV.File("data/growth.csv"))

avg_df = DataFrame(Time = df.Time)
avg_df.gluc = (df.C4 + df.D4 + df.E4) / 3
avg_df.no_cys = (df.C5 + df.D5 + df.E5) / 3
avg_df.no_pro = (df.C6 + df.D6 + df.E6) / 3
avg_df.no_arg = (df.C7 + df.D7 + df.E7) / 3
avg_df.no_leu = (df.C8 + df.D8 + df.E8) / 3
avg_df.no_val = (df.C9 + df.D9 + df.E9) / 3
avg_df.ribose = (df.C10 + df.D10 + df.E10) / 3
xs = 1:97

f = Figure(size = (800, 400))
ax = Axis(
    f[1,1];
    ylabel = "OD600",
    xlabel = "Time point"
)
lines!(ax, xs, avg_df.gluc, label = "Glucose + 5AAs", color=:blue)
lines!(ax, xs, avg_df.no_cys, label = "no Cys", color=:red)
lines!(ax, xs, avg_df.no_pro, label = "no Pro", color=:green)
lines!(ax, xs, avg_df.no_arg, label = "no Arg", color=:purple)
lines!(ax, xs, avg_df.no_leu, label = "no Leu", color=:orange)
lines!(ax, xs, avg_df.no_val, label = "no Val", color=:black)
lines!(ax, xs, avg_df.ribose, label = "Ribose + 5AAs",color=:brown)
f
f[1, 2] = Legend(f, ax, "Media", framevisible = false)

f

save("data/growth_curve.png",f)
