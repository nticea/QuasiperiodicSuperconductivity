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
files = readdir(joinpath(@__DIR__, "data"))
df = load_dfs()
df = df[(df.T.>=T_cutoff), :]

Js = sort(unique(df.J))
Ls = sort(unique(df.L))
cmap = cgrad(:matter, length(Ls), categorical=true)
mmap = [:dashdot, :dash, :solid]

gdf = groupby(df[(df.J.==(Js[2])), :], [:L])
ps, pd = plot(), plot()
psχ, pdχ = plot(), plot()
for g in gdf
    # get the data 
    L, J = g.L[1], g.J[1]
    @show L, J
    Lidx = findall(x -> x == L, Ls)[1]
    Jidx = findall(x -> x == J, Js)[1]

    Ts, dχs, χs = g.T, g.dχdlogT, g.χ
    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    dχs = dχs[sortidx]
    χs = χs[sortidx]
    dχswave, dχdwave = [], []
    χswave, χdwave = [], []
    for (χ, dχ) in zip(χs, dχs)
        s, d = uniform_susceptibility_components(χ, ndims=ndims)
        push!(χswave, s)
        push!(χdwave, d)

        s, d = uniform_susceptibility_components(dχ, ndims=ndims)
        push!(dχswave, s)
        push!(dχdwave, d)
    end

    ps = plot!(ps, Ts, dχswave, c=cmap[Lidx], m=mmap[Jidx], label=nothing, xaxis=:log10)
    ps = scatter!(ps, Ts, dχswave, c=cmap[Lidx], m=mmap[Jidx], label="J=$J, L=$L", xaxis=:log10)
    pd = plot!(pd, Ts, dχdwave, c=cmap[Lidx], m=mmap[Jidx], label=nothing, xaxis=:log10)
    pd = scatter!(pd, Ts, dχdwave, c=cmap[Lidx], m=mmap[Jidx], label="J=$J, L=$L", xaxis=:log10)

    psχ = plot!(psχ, Ts, χswave, c=cmap[Lidx], m=mmap[Jidx], label=nothing, xaxis=:log10)
    psχ = scatter!(psχ, Ts, χswave, c=cmap[Lidx], m=mmap[Jidx], label="J=$J, L=$L", xaxis=:log10)
    pdχ = plot!(pdχ, Ts, χdwave, c=cmap[Lidx], m=mmap[Jidx], label=nothing, xaxis=:log10)
    pdχ = scatter!(pdχ, Ts, χdwave, c=cmap[Lidx], m=mmap[Jidx], label="J=$J, L=$L", xaxis=:log10)
end