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

ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
t = 1
μ = 1e-8
ϕx, ϕy, ϕz = 0, 0, 0
savefigs = false
T_cutoff = 1e-2

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
df = load_dfs()
Ls = unique(df.L)

if length(Ls) == 1
    cmap = ["red"]
else
    cmap = cgrad(:matter, length(Ls), categorical=true)
end

px, pos = plot(), plot()
for (l, L) in enumerate(Ls)
    # extract only the parameters we are interested in 
    dfL = df[(df.L.==L).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy).&(df.ϕz.==ϕz).&(df.θ.==θ).&(df.Q.==Q).&(df.ndims.==ndims), :]
    Js = sort(unique(dfL.J))

    for (j, J) in enumerate(Js)
        dfJ = dfL[(dfL.J.==J), :]
        Ts = dfJ.T
        χswave = dfJ.χswave
        χdwave = dfJ.χdwave

        sortidx = sortperm(Ts)
        Ts = Ts[sortidx]
        χswave = χswave[sortidx]
        χdwave = χdwave[sortidx]

        # consider only the linear regime 
        # χs = χs[Ts.>=T_cutoff]
        # Ts = Ts[Ts.>=T_cutoff]

        if length(Ts) > 0
            #plot!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10)
            scatter!(px, Ts, χdwave, color=cmap[l], label="J=$J, L=$L", xaxis=:log10)

            #plot!(pos, Ts, χswave, color=cmap[j], label=nothing, xaxis=:log10)
            scatter!(pos, Ts, χswave, color=cmap[l], label="J=$J, L=$L", xaxis=:log10)

            title!(px, "Susceptibility (d-wave)")
            xlabel!(px, "T")
            ylabel!(px, "χ")
            title!(pos, "Susceptibility (s-wave)")
            xlabel!(pos, "T")
            ylabel!(pos, "χ")
        end
    end
end

p = plot(px, pos, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility")
# if savefigs
#     savefig(p, joinpath(figpath, "susceptibility.pdf"))
# end