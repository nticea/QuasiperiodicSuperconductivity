## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using CSV
using DataFrames
using StatsPlots

include("../../src/model.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")

## PARAMETERS ## 

L = 23 # the full system is L × L 
Q = (√5 - 1) / 2
θ = π / 7
savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
files = readdir(joinpath(@__DIR__, "data"))
df = DataFrame(L=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[],
    T=[], χ=[])
for f in files
    if endswith(f, ".csv")
        dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
        append!(df, dfi)
    end
end

# extract only the parameters we are interested in 
df = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q), :]
Js = sort(unique(df.J))
Js = Js[1:2:end]

ps = []
cmap = cgrad(:matter, length(Js), categorical=true)
px, py, pos = plot(), plot(), plot()
for (j, J) in enumerate(Js)
    dfJ = df[(df.J.==J), :]
    Ts = dfJ.T
    χs = dfJ.χ
    if length(Ts) > 0

        # on-site
        χswave = 1 / (L * L) * [a[9] for a in χs]

        # make the d-wave components 
        xx = [a[1] for a in χs]
        yy = [a[5] for a in χs]
        xy = [a[2] for a in χs]
        yx = [a[4] for a in χs]
        χdwave = 1 / (L * L) * (xx + yy - xy - yx)

        plot!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(px, Ts, χdwave, color=cmap[j], label="J=$J", xaxis=:log10)
        plot!(pos, Ts, χswave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(pos, Ts, χswave, color=cmap[j], label="J=$J", xaxis=:log10)

        title!(px, "Susceptibility (d-wave)")
        xlabel!(px, "T")
        ylabel!(px, "χ")
        title!(pos, "Susceptibility (s-wave)")
        xlabel!(pos, "T")
        ylabel!(pos, "χ")
    end
end

p = plot(ps..., layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="LGE χ for $L × $L lattice")
if savefigs
    savefig(p, joinpath(figpath, "susceptibility_$(L)L.pdf"))
end