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
J_cutoff = 0.5

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
df = load_dfs(dirname="data_PBC")
df = df[(df.T.>=T_cutoff).&&(df.J.<=J_cutoff), :]

Js = sort(unique(df.J))
Ls = sort(unique(df.L))
cmap = cgrad(:matter, length(Js), categorical=true)
mmap = [:solid, :dash]

ps, pd = plot(title="swave disorder potential at μ=0.5", xlabel="T", ylabel="dχ/dlogT"), plot(title="dwave disorder potential at μ=0.5", xlabel="T", ylabel="dχ/dlogT", legend=:right)
psχ, pdχ = plot(title="swave disorder potential at μ=0.5", xlabel="T", ylabel="χ"), plot(title="dwave disorder potential at μ=0.5", xlabel="T", ylabel="χ")
for (Jidx, J) in enumerate(Js)
    for (Lidx, L) in enumerate(Ls)
        # get the data 
        g = df[(df.J.==J).&&(df.L.==L), :]

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

        ps = plot!(ps, Ts, dχswave, c=cmap[Jidx], ls=mmap[Lidx], label="J=$J, L=$L", xaxis=:log10)
        ps = scatter!(ps, Ts, dχswave, c=cmap[Jidx], label=nothing, xaxis=:log10)
        pd = plot!(pd, Ts, dχdwave, c=cmap[Jidx], ls=mmap[Lidx], label="J=$J, L=$L", xaxis=:log10)
        pd = scatter!(pd, Ts, dχdwave, c=cmap[Jidx], label=nothing, xaxis=:log10)

        psχ = plot!(psχ, Ts, χswave, c=cmap[Jidx], ls=mmap[Lidx], label="J=$J, L=$L", xaxis=:log10)
        psχ = scatter!(psχ, Ts, χswave, c=cmap[Jidx], label=nothing, xaxis=:log10)
        pdχ = plot!(pdχ, Ts, χdwave, c=cmap[Jidx], ls=mmap[Lidx], label="J=$J, L=$L", xaxis=:log10)
        pdχ = scatter!(pdχ, Ts, χdwave, c=cmap[Jidx], label=nothing, xaxis=:log10)
    end
end

