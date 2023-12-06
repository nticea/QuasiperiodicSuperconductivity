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
    dirname = "3D_data_random_PBC"
    title = "disordered potential"
else
    dirname = "3D_data_PBC"
    title = "quasiperiodic potential"
end
df = load_dfs(dirname=dirname)
df = df[(df.L.==L).&&(df.J.>J_cutoff).&&(df.T.<=T_max).&&(df.T.>=T_min).&&(df.Q.==Q).&&(df.θ.==θ).&&(df.ndims.==ndims), :]

T_max_act = maximum(df.T)
T_min_act = minimum(df.T)

# make a new dataframe to hold the averaged quantities 
gdf = groupby(df, [:J, :T])
dfsummary = DataFrame(J=[], T=[], χswave=[], χdwave=[], dχswave=[], dχdwave=[])
for g in gdf
    J, T = g.J[1], g.T[1]
    χs, dχs = g.χ, g.dχdlogT
    for (χ, dχ) in zip(χs, dχs)
        if ndims == 2
            χ = reshape(χ, (3, 3))
            dχ = reshape(dχ, (3, 3))
        elseif ndims == 3
            χ = reshape(χ, (4, 4))
            dχ = reshape(dχ, (4, 4))
        else
            println("sucks to suck")
        end
        χsw, χdw = uniform_susceptibility_components(χ, ndims=ndims)
        dχsw, dχdw = uniform_susceptibility_components(dχ, ndims=ndims)
        dfi = DataFrame(J=[J], T=[T], χswave=[χsw], χdwave=[χdw], dχswave=[dχsw], dχdwave=[dχdw])
        append!(dfsummary, dfi)
    end
end
gdf = groupby(dfsummary, [:J, :T])
dfmean = combine(gdf, [:χswave => mean, :χdwave => mean, :dχswave => mean, :dχdwave => mean,
    :χswave => sem, :χdwave => sem, :dχswave => sem, :dχdwave => sem])

sems = dfmean.dχdwave_sem
sems[isnan.(sems)] .= 0
dfmean.dχdwave_sem = sems

sems = dfmean.dχswave_sem
sems[isnan.(sems)] .= 0
dfmean.dχswave_sem = sems

Js = sort(unique(dfmean.J))
Ts = sort(unique(dfmean.T))
cmap = reverse(cgrad(:viridis, length(Ts), categorical=true))
#cmap = ["red", "blue", "green", "orange"]
pχswave, pχdwave, pdχswave, pdχdwave = plot(title="χ swave for $title", grid=false), plot(title="χ dwave for $title", grid=false), plot(title="dχdlogT swave for $title", grid=false), plot(title="dχdlogT dwave for $title", grid=false)
for (Tᵢ, T) in enumerate(reverse(Ts))
    dfT = dfmean[(dfmean.T.==T), :]
    Js, χswave, χdwave, dχswave, dχdwave = dfT.J, dfT.χswave_mean, dfT.χdwave_mean, dfT.dχswave_mean, dfT.dχdwave_mean
    σ_χswave, σ_χdwave, σ_dχswave, σ_dχdwave = dfT.χswave_sem, dfT.χdwave_sem, dfT.dχswave_sem, dfT.dχdwave_sem
    sortidx = sortperm(Js)
    Js = Js[sortidx]
    χswave = χswave[sortidx]
    χdwave = χdwave[sortidx]
    dχswave = dχswave[sortidx]
    dχdwave = dχdwave[sortidx]
    σ_χswave = σ_χswave[sortidx]
    σ_χdwave = σ_χdwave[sortidx]
    σ_dχswave = σ_dχswave[sortidx]
    σ_dχdwave = σ_dχdwave[sortidx]

    plot!(pχswave, Js, χswave, xlabel="J", ylabel="χ", label=nothing, c=cmap[Tᵢ])
    scatter!(pχswave, Js, χswave, xlabel="J", ylabel="χ", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pχswave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pχdwave, Js, χdwave, xlabel="J", ylabel="χ", label=nothing, c=cmap[Tᵢ])
    scatter!(pχdwave, Js, χdwave, xlabel="J", ylabel="χ", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pχdwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pdχswave, Js, dχswave, xlabel="J", ylabel="dχ/dT", label=nothing, c=cmap[Tᵢ])
    scatter!(pdχswave, Js, dχswave, xlabel="J", ylabel="dχ/dT", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pdχswave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pdχdwave, Js, dχdwave, xlabel="J", ylabel="dχ/dT", label=nothing, c=cmap[Tᵢ])
    scatter!(pdχdwave, Js, dχdwave, xlabel="J", ylabel="dχ/dT", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pdχdwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)
end

plot!(pdχdwave, xlims=[0, 1.5], ylims=[-0.2, 0], legend=false)
plot!(pdχswave, xlims=[0, 1.5], ylims=[-0.2, 0], legend=false)
plot!(pχdwave, xlims=[0, 1.5], legend=false)
plot!(pχswave, xlims=[0, 1.5], legend=false)