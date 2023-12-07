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

L = 11 # the full system is L × L 
μ = 0.75
ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
T_max = 0.1
T_min = 0#0.002
disorder = false
J_cutoff = -1

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
if disorder
    dirname = "data_random_3D"
    title = "disordered potential"
else
    dirname = "data_QP_3D"
    title = "quasiperiodic potential"
end
df = load_dfs(dirname=dirname)
c = df[(df.L.==L).&&(df.J.>J_cutoff).&&(df.T.<=T_max).&&(df.T.>=T_min).&&(df.Q.==Q).&&(df.θ.==θ).&&(df.ndims.==ndims), :]

T_max_act = maximum(df.T)
T_min_act = minimum(df.T)

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

# sems = dfmean.dχdwave_sem
# sems[isnan.(sems)] .= 0
# dfmean.dχdwave_sem = sems

# sems = dfmean.dχswave_sem
# sems[isnan.(sems)] .= 0
# dfmean.dχswave_sem = sems

Js = sort(unique(dfmean.J))
Ts = sort(unique(dfmean.T))
cmap = reverse(cgrad(:viridis, length(Ts), categorical=true))
#cmap = ["red", "blue", "green", "orange"]
pχswave, pdχswave = plot(title="χ swave for $title", grid=false), plot(title="dχdlogT swave for $title", grid=false)
pχdwave, pdχdwave = plot(title="χ dwave for $title", grid=false), plot(title="dχdlogT dwave for $title", grid=false)
pχpwave, pdχpwave = plot(title="χ pwave for $title", grid=false), plot(title="dχdlogT pwave for $title", grid=false)

for (Tᵢ, T) in enumerate(reverse(Ts))
    dfT = dfmean[(dfmean.T.==T), :]
    Js, χswave, χdwave, χpwave, dχswave, dχdwave, dχpwave = dfT.J, dfT.χswave_mean, dfT.χdwave_mean, dfT.χpwave_mean, dfT.dχswave_mean, dfT.dχdwave_mean, dfT.dχpwave_mean
    σ_χswave, σ_χdwave, σ_χpwave, σ_dχswave, σ_dχdwave, σ_dχpwave = dfT.χswave_sem, dfT.χdwave_sem, dfT.χpwave_sem, dfT.dχswave_sem, dfT.dχdwave_sem, dfT.dχpwave_sem
    sortidx = sortperm(Js)
    Js = Js[sortidx]
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

    plot!(pχswave, Js, χswave, xlabel="J", ylabel="χ", label=nothing, c=cmap[Tᵢ])#, ribbon=σ_χswave)
    scatter!(pχswave, Js, χswave, xlabel="J", ylabel="χ", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pχswave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pχdwave, Js, χdwave, xlabel="J", ylabel="χ", label=nothing, c=cmap[Tᵢ])#, ribbon=σ_χdwave)
    scatter!(pχdwave, Js, χdwave, xlabel="J", ylabel="χ", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pχdwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pχpwave, Js, χpwave, xlabel="J", ylabel="χ", label=nothing, c=cmap[Tᵢ])#, ribbon=σ_χpwave)
    scatter!(pχpwave, Js, χpwave, xlabel="J", ylabel="χ", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pχpwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pdχswave, Js, dχswave, xlabel="J", ylabel="dχ/dT", label=nothing, c=cmap[Tᵢ])#, ribbon=σ_dχswave)
    scatter!(pdχswave, Js, dχswave, xlabel="J", ylabel="dχ/dT", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pdχswave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pdχdwave, Js, dχdwave, xlabel="J", ylabel="dχ/dT", label=nothing, c=cmap[Tᵢ])#, ribbon=σ_dχdwave)
    scatter!(pdχdwave, Js, dχdwave, xlabel="J", ylabel="dχ/dT", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pdχdwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)

    plot!(pdχpwave, Js, dχpwave, xlabel="J", ylabel="dχ/dT", label=nothing, c=cmap[Tᵢ])#, ribbon=σ_dχpwave)
    scatter!(pdχpwave, Js, dχpwave, xlabel="J", ylabel="dχ/dT", label="T=$T", c=cmap[Tᵢ])
    heatmap!(pdχpwave, [T_min_act T_min_act; T_max_act T_max_act], cmap=:viridis, clims=(T_min_act, T_max_act), alpha=0)
end

plot!(pdχdwave, xlims=[0, 1.5], ylims=[-0.2, 0], legend=false)
plot!(pdχswave, xlims=[0, 1.5], ylims=[-0.2, 0], legend=false)
plot!(pdχpwave, xlims=[0, 1.5], ylims=[-0.2, 0], legend=false)
plot!(pχdwave, xlims=[0, 1.5], legend=false)
plot!(pχswave, xlims=[0, 1.5], legend=false)
plot!(pχpwave, xlims=[0, 1.5], legend=false)