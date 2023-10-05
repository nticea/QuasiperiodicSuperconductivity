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

L = 17 # the full system is L × L 
ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
T_cutoff = 0#1e-1

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
dirname = "data_PBC"
df = load_dfs(dirname=dirname)
df = df[(df.L.==L).&&(df.T.>=T_cutoff).&&(df.Q.==Q).&&(df.θ.==θ).&&(df.ndims.==ndims), :]

# make a new dataframe to hold the averaged quantities 
gdf = groupby(df, [:J, :T])
dfsummary = DataFrame(J=[], T=[], χswave=[], χdwave=[], dχswave=[], dχdwave=[])
function mean_χ(χs)
    χs = hcat(χs...)
    return mean(χs, dims=2)
end
for g in gdf
    J, T = g.J[1], g.T[1]
    χs, dχs = g.χ, g.dχdlogT
    for (χ, dχ) in zip(χs, dχs)
        χ = reshape(χ, (4, 4))
        dχ = reshape(dχ, (4, 4))
        χsw, χdw = uniform_susceptibility_components(χ, ndims=ndims)
        dχsw, dχdw = uniform_susceptibility_components(dχ, ndims=ndims)
        dfi = DataFrame(J=[J], T=[T], χswave=[χsw], χdwave=[χdw], dχswave=[dχsw], dχdwave=[dχdw])
        append!(dfsummary, dfi)
    end
end
gdf = groupby(dfsummary, [:J, :T])
dfmean = combine(gdf, [:χswave => mean, :χdwave => mean, :dχswave => mean, :dχdwave => mean])


Js = sort(unique(dfmean.J))
cmap = cgrad(:matter, length(Js), categorical=true)
pχswave, pχdwave, pdχswave, pdχdwave = plot(title="χ swave"), plot(title="χ dwave"), plot(title="dχdlogT swave"), plot(title="dχdlogT dwave")
for J in Js
    dfJ = dfmean[(dfmean.J.==J), :]
    Ts, χswave, χdwave, dχswave, dχdwave = dfmean.T, dfmean.χswave_mean, dfmean.χdwave_mean, dfmean.dχswave_mean, dfmean.dχdwave_mean
end




