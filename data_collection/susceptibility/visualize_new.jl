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

L = 23 # the full system is L × L 
μ = 0.75
ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
T_max = 0.1
T_min = 0.002
disorder = false
J_cutoff = -1

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
if disorder
    # dirname = "data_random_3D"
    dirname = "FINAL_3D_data_random_PBC"
    title = "disordered potential"
else
    dirname = "FINAL_3D_data_PBC"
    title = "quasiperiodic potential"
end
df = load_dfs(dirname=dirname)
c = df[(df.L.==L).&&(df.J.>J_cutoff).&&(df.T.<=T_max).&&(df.T.>=T_min).&&(df.Q.==Q).&&(df.θ.==θ).&&(df.ndims.==ndims), :]

T_max_act = maximum(df.J)
T_min_act = minimum(df.J)

# make a new dataframe to hold the averaged quantities 
gdf = groupby(df, [:J, :T])
dfsummary = DataFrame(J=[], T=[], χswave=[], χdwave=[], χpwave=[], dχswave=[], dχdwave=[], dχpwave=[])
for g in gdf
    J, T = g.J[1], g.T[1]
    χs, dχs = g.χ_dwave, g.dχdlogT_dwave
    χps, dχps = g.χ_pwave, g.dχdlogT_pwave
    for (χ, dχ, χp, dχp) in zip(χs, dχs, χps, dχps)
        if ndims == 2
            χ = reshape(χ, (3, 3))
            dχ = reshape(dχ, (3, 3))
            χp = reshape(χp, (3, 3))
            dχp = reshape(dχp, (3, 3))
        elseif ndims == 3
            χ = reshape(χ, (4, 4))
            dχ = reshape(dχ, (4, 4))
            χp = reshape(χp, (4, 4))
            dχp = reshape(dχp, (4, 4))
        else
            println("sucks to suck")
        end
        χsw, χdw = uniform_susceptibility_components(χ, ndims=ndims)
        dχsw, dχdw = uniform_susceptibility_components(dχ, ndims=ndims)
        _, χpw = uniform_susceptibility_components(χp, ndims=ndims)
        _, dχpw = uniform_susceptibility_components(dχp, ndims=ndims)
        dfi = DataFrame(J=[J], T=[T], χswave=[χsw], χdwave=[χdw], χpwave=[χpw], dχswave=[dχsw], dχdwave=[dχdw], dχpwave=[dχpw])
        append!(dfsummary, dfi)
    end
end
gdf = groupby(dfsummary, [:J, :T])
dfmean = combine(gdf, [:χswave => mean, :χdwave => mean, :χpwave => mean,
    :dχswave => mean, :dχdwave => mean, :dχpwave => mean,
    :χswave => sem, :χdwave => sem, :χpwave => sem,
    :dχswave => sem, :dχdwave => sem, :dχpwave => sem])

Js = sort(unique(dfmean.J))
cmap = reverse(cgrad(:viridis, length(Js), categorical=true))
#cmap = ["red", "blue", "green", "orange"]
pχswave, pdχswave = plot(title="χ swave for $title", grid=false), plot(title="dχdlogT swave for $title", grid=false)
pχdwave, pdχdwave = plot(title="χ dwave for $title", grid=false), plot(title="dχdlogT dwave for $title", grid=false)
pχpwave, pdχpwave = plot(title="χ pwave for $title", grid=false), plot(title="dχdlogT pwave for $title", grid=false)

for (Jᵢ, J) in enumerate(Js)
    dfJ = dfmean[(dfmean.J.==J), :]
    Ts, χswave, χdwave, χpwave, dχswave, dχdwave, dχpwave = dfJ.T, dfJ.χswave_mean, dfJ.χdwave_mean, dfJ.χpwave_mean, dfJ.dχswave_mean, dfJ.dχdwave_mean, dfJ.dχpwave_mean
    σ_χswave, σ_χdwave, σ_χpwave, σ_dχswave, σ_dχdwave, σ_dχpwave = dfJ.χswave_sem, dfJ.χdwave_sem, dfJ.χpwave_sem, dfJ.dχswave_sem, dfJ.dχdwave_sem, dfJ.dχpwave_sem
    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χswave = χswave[sortidx]
    χdwave = χdwave[sortidx]
    χpwave = χpwave[sortidx]
    dχswave = dχswave[sortidx]
    dχdwave = dχdwave[sortidx]
    dχpwave = dχpwave[sortidx]
    σ_χswave = σ_χswave[sortidx]
    σ_χdwave = σ_χdwave[sortidx]
    σ_χpwave = σ_χpwave[sortidx]
    σ_dχswave = σ_dχswave[sortidx]
    σ_dχdwave = σ_dχdwave[sortidx]
    σ_dχpwave = σ_dχpwave[sortidx]

    plot!(pχswave, Ts, χswave, xlabel="T", ylabel="χ", label=nothing, c=cmap[Jᵢ], ribbon=σ_χswave)
    scatter!(pχswave, Ts, χswave, xlabel="T", ylabel="χ", label="J=$J", c=cmap[Jᵢ])
    #heatmap!(pχswave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pχdwave, Ts, χdwave, xlabel="T", ylabel="χ", label=nothing, c=cmap[Jᵢ], ribbon=σ_χdwave)
    scatter!(pχdwave, Ts, χdwave, xlabel="T", ylabel="χ", label="J=$J", c=cmap[Jᵢ])
    #heatmap!(pχdwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pχpwave, Ts, χpwave, xlabel="T", ylabel="χ", label=nothing, c=cmap[Jᵢ], ribbon=σ_χpwave)
    scatter!(pχpwave, Ts, χpwave, xlabel="T", ylabel="χ", label="J=$J", c=cmap[Jᵢ])
    #heatmap!(pχpwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pdχswave, Ts, dχswave, xlabel="T", ylabel="dχ/dT", label=nothing, c=cmap[Jᵢ], ribbon=σ_dχswave)
    scatter!(pdχswave, Ts, dχswave, xlabel="T", ylabel="dχ/dT", label="J=$J", c=cmap[Jᵢ])
    #heatmap!(pdχswave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pdχdwave, Ts, dχdwave, xlabel="T", ylabel="dχ/dT", label=nothing, c=cmap[Jᵢ], ribbon=σ_dχdwave)
    scatter!(pdχdwave, Ts, dχdwave, xlabel="T", ylabel="dχ/dT", label="J=$J", c=cmap[Jᵢ])
    #heatmap!(pdχdwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pdχpwave, Ts, dχpwave, xlabel="T", ylabel="dχ/dT", label=nothing, c=cmap[Jᵢ], ribbon=σ_dχpwave)
    scatter!(pdχpwave, Ts, dχpwave, xlabel="T", ylabel="dχ/dT", label="J=$J", c=cmap[Jᵢ])
    #heatmap!(pdχpwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)
end

plot!(pdχdwave, ylims=[-0.2, 0], legend=false, xaxis=:log10)
plot!(pdχswave, ylims=[-0.2, 0], legend=false, xaxis=:log10)
plot!(pdχpwave, ylims=[-0.2, 0], legend=false, xaxis=:log10)
plot!(pχdwave, legend=false, xaxis=:log10)
plot!(pχswave, legend=false, xaxis=:log10)
plot!(pχpwave, legend=false, xaxis=:log10)