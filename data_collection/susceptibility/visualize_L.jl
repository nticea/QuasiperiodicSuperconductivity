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
savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
T_cutoff = 1e-2

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
files = readdir(joinpath(@__DIR__, "data"))
df = load_dfs()
Js = sort(unique(df.J))
Ls = sort(unique(df.L))

if length(Ls) == 1
    cmap = ["red"]
else
    cmap = cgrad(:matter, length(Ls), categorical=true)
end

ps = []
for J in Js
    px, pos = plot(), plot()
    for (l, L) in enumerate(Ls)
        # extract only the parameters we are interested in 
        dfL = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q), :]

        dfJ = dfL[(dfL.J.==J), :]
        Ts = dfJ.T
        χs = dfJ.χ
        χs = [reshape(χ, 4, 4) for χ in χs]

        sortidx = sortperm(Ts)
        Ts = Ts[sortidx]
        χs = χs[sortidx]

        if length(Ts) > 0

            # on-site
            χswave = [χ[1, 1] for χ in χs]

            # make the d-wave components 
            xx, yy = [χ[2, 2] for χ in χs], [χ[3, 3] for χ in χs]
            xy, yx = [χ[2, 3] for χ in χs], [χ[3, 2] for χ in χs]
            if ndims == 2
                χdwave = xx + yy - xy - yx
            elseif ndims == 3
                zz = [χ[4, 4] for χ in χs]
                xz, zx = [χ[2, 4] for χ in χs], [χ[4, 2] for χ in χs]
                yz, zy = [χ[3, 4] for χ in χs], [χ[4, 3] for χ in χs]
                χdwave = xx + yy + zz - xy - yx - xz - zx - yz - zy
            else
                println("sorry")
                χdwave = nothing
            end

            plot!(px, Ts, χdwave, color=cmap[l], label=nothing, xaxis=:log10)
            scatter!(px, Ts, χdwave, color=cmap[l], label="J=$J, L=$L", xaxis=:log10)

            plot!(pos, Ts, χswave, color=cmap[l], label=nothing, xaxis=:log10)
            scatter!(pos, Ts, χswave, color=cmap[l], label="J=$J, L=$L", xaxis=:log10)

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
    push!(ps, p)
end
# if savefigs
#     savefig(p, joinpath(figpath, "susceptibility.pdf"))
# end