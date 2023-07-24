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

L = 11 # the full system is L × L 
Q = (√5 - 1) / 2
θ = π / 7
savefigs = false

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
files = readdir(joinpath(@__DIR__, "data"))
df = DataFrame(L=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[],
    symmetry=[], T=[], χ=[])
for f in files
    if endswith(f, ".csv")
        dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
        # process all of the arrays 
        if dfi.symmetry[1] == "d-wave"
            dfi = convert_df_arrays(dfi, "χ")
        end
        append!(df, dfi)
    end
end

# extract only the parameters we are interested in 
df = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q), :]
Js = sort(unique(df.J))

ps = []
for symmetry in ["d-wave"]
    cmap = cgrad(:matter, length(Js), categorical=true)
    px, py, pos = plot(), plot(), plot()
    for (j, J) in enumerate(Js)
        dfJ = df[(df.symmetry.==symmetry).&(df.J.==J), :]
        Ts = dfJ.T
        χs = dfJ.χ
        if length(Ts) > 0

            χx = [a[1] for a in χs]
            χy = [a[5] for a in χs]
            χos = [a[9] for a in χs]
            plot!(px, Ts, χx, color=cmap[j], label=nothing, xaxis=:log10)
            scatter!(px, Ts, χx, color=cmap[j], label="J=$J", xaxis=:log10)
            plot!(py, Ts, χy, color=cmap[j], label=nothing, xaxis=:log10)
            scatter!(py, Ts, χy, color=cmap[j], label="J=$J", xaxis=:log10)
            plot!(pos, Ts, χos, color=cmap[j], label=nothing, xaxis=:log10)
            scatter!(pos, Ts, χos, color=cmap[j], label="J=$J", xaxis=:log10)

            title!(px, "Susceptibility (x̂)")
            xlabel!(px, "T")
            ylabel!(px, "χ")
            title!(py, "Susceptibility (ŷ)")
            xlabel!(py, "T")
            ylabel!(py, "χ")
            title!(pos, "Susceptibility (on-site)")
            xlabel!(pos, "T")
            ylabel!(pos, "χ")
        end
    end
    push!(ps, px, py, pos)
end

p = plot(ps..., layout=Plots.grid(1, 3,
        widths=[1 / 3, 1 / 3, 1 / 3]), size=(1700, 500), plot_title="LGE χ for $L × $L lattice")
