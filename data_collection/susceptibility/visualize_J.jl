## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots, StatsPlots
using CSV, DataFrames
using CurveFit

include("../../src/model.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")
include("utilities.jl")

## PARAMETERS ## 

L = 105 # the full system is L × L 
Q = (√5 - 1) / 2
θ = π / 7
savefigs = true
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
        dfi = convert_df_arrays(dfi, "χ")
        append!(df, dfi)
    end
end

px, pos = plot(), plot()
# extract only the parameters we are interested in 
dfL = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q), :]
Js = sort(unique(dfL.J))
cmap = cgrad(:matter, length(Js), categorical=true)

for (j, J) in enumerate(Js)
    dfJ = dfL[(dfL.J.==J), :]
    Ts = dfJ.T
    χs = dfJ.χ

    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χs = χs[sortidx]

    if length(Ts) > 0

        # on-site
        χswave = [a[9] for a in χs]

        # make the d-wave components 
        xx = [a[1] for a in χs]
        yy = [a[5] for a in χs]
        xy = [a[2] for a in χs]
        yx = [a[4] for a in χs]
        χdwave = 1 / 4 * (xx + yy - xy - yx)

        plot!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(px, Ts, χdwave, color=cmap[j], label="J=$J, L=$L", xaxis=:log10)

        plot!(pos, Ts, χswave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(pos, Ts, χswave, color=cmap[j], label="J=$J, L=$L", xaxis=:log10)

        title!(px, "Susceptibility (d-wave)")
        xlabel!(px, "T")
        ylabel!(px, "χ")
        title!(pos, "Susceptibility (s-wave)")
        xlabel!(pos, "T")
        ylabel!(pos, "χ")
    end
end

p = plot(px, pos, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility")
if savefigs
    savefig(p, joinpath(figpath, "susceptibility.pdf"))
end