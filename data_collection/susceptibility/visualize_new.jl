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
T_cutoff = 0#2 * 1e-3
disorder = true

# read files 
if savefigs
    figpath = mkpath(joinpath(@__DIR__, "figures"))
end
if disorder
    dirname = "3D_data_random_PBC"
    title = "disordered potential"
else
    dirname = "3D_data_PBC"
    title = "quasiperiodic potential"
end
# if disorder
#     dirname = "data_random_PBC"
#     title = "disordered potential"
#     pot = "disorder"
# else
#     dirname = "data_PBC"
#     title = "quasiperiodic potential"
#     pot = "QP"
# end
df = load_dfs(dirname=dirname)
df = df[(df.L.==L).&&(df.T.>=T_cutoff).&&(df.Q.==Q).&&(df.θ.==θ).&&(df.ndims.==ndims), :]

# xtcks = LinRange(T_cutoff, 10, 20)
xtcks = [1e-3, 1e-2, 1e-1, 1, 10]
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
cmap = cgrad(:matter, length(Js), categorical=true)
#cmap = ["red", "blue", "green", "orange"]
pχswave, pχdwave, pdχswave, pdχdwave = plot(title="χ swave for $title", grid=false, ylims=(0, 1.5), xticks=(xtcks, xtcks)), plot(title="χ dwave for $title", grid=false, ylims=(0, 1.5), xticks=(xtcks, xtcks)), plot(title="dχdlogT swave for $title", grid=false, ylims=(-0.25, 0), xticks=(xtcks, xtcks)), plot(title="dχdlogT dwave for $title", grid=false, ylims=(-0.25, 0), xticks=(xtcks, xtcks))
for (Jᵢ, J) in enumerate(Js)
    dfJ = dfmean[(dfmean.J.==J), :]
    Ts, χswave, χdwave, dχswave, dχdwave = dfJ.T, dfJ.χswave_mean, dfJ.χdwave_mean, dfJ.dχswave_mean, dfJ.dχdwave_mean
    σ_χswave, σ_χdwave, σ_dχswave, σ_dχdwave = dfJ.χswave_sem, dfJ.χdwave_sem, dfJ.dχswave_sem, dfJ.dχdwave_sem
    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χswave = χswave[sortidx]
    χdwave = χdwave[sortidx]
    dχswave = dχswave[sortidx]
    dχdwave = dχdwave[sortidx]
    σ_χswave = σ_χswave[sortidx]
    σ_χdwave = σ_χdwave[sortidx]
    σ_dχswave = σ_dχswave[sortidx]
    σ_dχdwave = σ_dχdwave[sortidx]

    plot!(pχswave, Ts, χswave, xaxis=:log10, xlabel="T", ylabel="χ", label=nothing, c=cmap[Jᵢ], ribbon=σ_χswave)
    scatter!(pχswave, Ts, χswave, xaxis=:log10, xlabel="T", ylabel="χ", label="J=$J", c=cmap[Jᵢ])

    plot!(pχdwave, Ts, χdwave, xaxis=:log10, xlabel="T", ylabel="χ", label=nothing, c=cmap[Jᵢ], ribbon=σ_χdwave)
    scatter!(pχdwave, Ts, χdwave, xaxis=:log10, xlabel="T", ylabel="χ", label="J=$J", c=cmap[Jᵢ])

    plot!(pdχswave, Ts, dχswave, xaxis=:log10, xlabel="T", ylabel="dχ", label=nothing, c=cmap[Jᵢ], ribbon=σ_dχswave)
    scatter!(pdχswave, Ts, dχswave, xaxis=:log10, xlabel="T", ylabel="dχ", label="J=$J", c=cmap[Jᵢ])

    plot!(pdχdwave, Ts, dχdwave, xaxis=:log10, xlabel="T", ylabel="dχ", label=nothing, c=cmap[Jᵢ], ribbon=σ_dχdwave)
    scatter!(pdχdwave, Ts, dχdwave, xaxis=:log10, xlabel="T", ylabel="dχ", label="J=$J", c=cmap[Jᵢ])
end

if savefigs
    savefig(pχswave, joinpath(figpath, "χswave_$(L)L_$(pot).pdf"))
    savefig(pχdwave, joinpath(figpath, "χdwave_$(L)L_$(pot).pdf"))
    savefig(pdχswave, joinpath(figpath, "dχswave_$(L)L_$(pot).pdf"))
    savefig(pdχdwave, joinpath(figpath, "dχdwave_$(L)L_$(pot).pdf"))
end
