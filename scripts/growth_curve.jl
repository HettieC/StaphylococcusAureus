using CSV, DataFrames, CairoMakie

df = DataFrame(CSV.File("data/experimental/hettie_2nd_trail.csv"))

avg_df = DataFrame(Time = df.Duration)
avg_df.gluc = (df.C4 + df.D4 + df.E4) / 3
avg_df.no_cys = (df.C5 + df.D5 + df.E5) / 3
avg_df.no_pro = (df.C6 + df.D6 + df.E6) / 3
avg_df.no_arg = (df.C7 + df.D7 + df.E7) / 3
avg_df.no_leu = (df.C8 + df.D8 + df.E8) / 3
avg_df.no_val = (df.C9 + df.D9 + df.E9) / 3
avg_df.ribose = (df.C10 + df.D10 + df.E10) / 3
xs = 1:97

inch = 96
pt = 4/3
cm = inch / 2.54

set_theme!(figure_padding=3)

f = Figure(; size=(10cm, 6cm))#, backgroundcolor=:transparent)

ax = Axis(
    f[1,1];
    backgroundcolor=:transparent,
    ylabel = "OD600",
    xlabel = "Time (h)",
    xlabelsize=6pt,
    ylabelsize=6pt,
    xticklabelsize=5pt,
    yticklabelsize=5pt,
    xticks = [24,36,48,60,72],
    ygridvisible=false,
    xgridvisible=false,

)
lines!(ax, avg_df.Time.+24, avg_df.gluc, label = "Glucose + 5AAs", color=:blue)
lines!(ax, avg_df.Time.+24, avg_df.no_cys, label = "no Cys", color=:red)
lines!(ax, avg_df.Time.+24, avg_df.no_pro, label = "no Pro", color=:green)
lines!(ax, avg_df.Time.+24, avg_df.no_arg, label = "no Arg", color=:purple)
lines!(ax, avg_df.Time.+24, avg_df.no_leu, label = "no Leu", color=:orange)
lines!(ax, avg_df.Time.+24, avg_df.no_val, label = "no Val", color=:black)
lines!(ax, avg_df.Time.+24, avg_df.ribose, label = "Ribose + 5AAs",color=:brown)
f
leg = Legend(
    f[1,2], 
    ax, 
    halign = :left,
    valign = :top,
    framevisible = false,
    labelsize = 5pt,
    tellheight=true,
    tellwidth=false,
    margin = (0,0,0,0),
)

colsize!(f.layout,1,Aspect(1,1.4))
f
xlims!(ax,(24,72))
f
save("data/growth_curve.png",f)
save("data/experimental/growth_curve.png", f, px_per_unit = 1200/inch)

