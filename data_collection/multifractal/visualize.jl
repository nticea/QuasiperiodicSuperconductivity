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
λ = 1 / 5
E₀ = 0

savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end

df = load_dfs()
# @assert 1 == 0
# df = df[(df.Q.==Q).&(df.θ.==θ).&(df.λ.==λ).&(df.ndims.==ndims).&(df.E₀.==E₀), :]

# for a fixed L and J, we want to average over all ϕx and ϕy 
gs = groupby(df, [:L, :J])
res = combine(gs, [:α₀ => mean, :ipr_real => mean, :ipr_k => mean])

# plot everything 
Ls = unique(res.L)
cmap = cgrad(:viridis, length(Ls), categorical=true)
p = plot(xlabel="J", ylabel="α₀")
for (Lᵢ, L) in enumerate(Ls)
    dfi = res[(res.L.==L), :]
    Js, α₀s = dfi.J, dfi.α₀_mean
    sortidx = sortperm(Js)
    Js, α₀s = Js[sortidx], α₀s[sortidx]
    p = plot!(Js, α₀s, label=nothing, c=cmap[Lᵢ])
    p = scatter!(p, Js, α₀s, label="L=$L", c=cmap[Lᵢ])
end
# ylims!(p, (2.5, 7.5))
# xlims!(p, (1.7, 2.5))

# plot everything 
Ls = unique(res.L)
cmap = cgrad(:viridis, length(Ls), categorical=true)
p2 = plot(xlabel="J", ylabel="IPR (real space)")
for (Lᵢ, L) in enumerate(Ls)
    dfi = res[(res.L.==L), :]
    p2 = plot!(dfi.J, dfi.ipr_real_mean, label=nothing, c=cmap[Lᵢ])
    p2 = scatter!(p2, dfi.J, dfi.ipr_real_mean, label="L=$L", c=cmap[Lᵢ])
end