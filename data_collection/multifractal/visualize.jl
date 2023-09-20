## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots, StatsPlots
using CSV, DataFrames, Statistics

include("../../src/model.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")
include("utilities.jl")

## PARAMETERS ## 

ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
ℓ = 3
E₀ = 1

savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end

df = load_dfs()
df = df[(df.Q.==Q).&(df.θ.==θ).&(df.ℓ.==ℓ).&(df.ndims.==ndims).&(df.E₀.==E₀), :]

# for a fixed L and J, we want to average over all ϕx and ϕy 
gs = groupby(df, [:L, :J])
res = DataFrame(L=[], J=[], α₀=[])
for g in gs
    append!(res, DataFrame(L=[g.L[1]], J=[g.J[1]], α₀=[mean(g.α₀)]))
end

# plot everything 
Ls = unique(res.L)
cmap = cgrad(:viridis, length(Ls), categorical=true)
p = plot(xlabel="J", ylabel="α₀")
for (Lᵢ, L) in enumerate(Ls)
    dfi = res[(res.L.==L), :]
    p = plot!(dfi.J, dfi.α₀, label=nothing, c=cmap[Lᵢ])
    p = scatter!(p, dfi.J, dfi.α₀, label="L=$L", c=cmap[Lᵢ])
end