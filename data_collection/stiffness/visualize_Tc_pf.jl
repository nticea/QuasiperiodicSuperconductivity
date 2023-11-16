## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using CSV
using DataFrames
using StatsPlots, Statistics

include("../../src/BdG.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")
include("utilities.jl")

## PARAMETERS ## 

L = 7 # the full system is L × L 
t = 1
Q = (√5 - 1) / 2
μ = 0.75
θ = π / 7
V0 = 1
V1 = -3
ndims = 3
periodic = true
disorder = true
savefigs = false
slice = 1

if disorder
    dirname = "data_$(ndims)D_disorder"
    pot = "disorder"
else
    dirname = "data_$(ndims)D_QP"
    pot = "QP"
end
if ndims == 3
    size_str = "$L × $L × $L"
elseif ndims == 2
    size_str = "$L × $L"
end

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
df_LGE_full = load_LGE(dirname)
df_BdG_full = load_BdG(dirname)

# for computing fermi velocities
m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=0, ϕy=0, ϕz=0, V0=0, V1=0, J=0, periodic=true, ndims=ndims, disorder=false)

# extract only the parameters we are interested in 
df_LGE = df_LGE_full[(df_LGE_full.L.==L).&(df_LGE_full.ndims.==ndims).&(df_LGE_full.θ.==θ).&(df_LGE_full.Q.==Q).&(df_LGE_full.V0.==V0).&(df_LGE_full.V1.==V1).&(df_LGE_full.μ.==μ), :]
df_BdG = df_BdG_full[(df_BdG_full.L.==L).&(df_BdG_full.ndims.==ndims).&(df_BdG_full.θ.==θ).&(df_BdG_full.Q.==Q).&(df_BdG_full.V0.==V0).&(df_BdG_full.V1.==V1).&(df_BdG_full.μ.==μ), :]

# Averaging
df_LGE = df_LGE[df_LGE.T.>0, :] # nonzero T
gdf_LGE = groupby(df_LGE, [:J])
LGE_mean = combine(gdf_LGE, [:λ => mean, :T => mean,
    :λ => sem, :T => sem])

# plotting
Js, Tcs, Tcs_err = LGE_mean.J, LGE_mean.T_mean, LGE_mean.T_sem
sortidx = sortperm(Js)
Js = Js[sortidx]
Tcs = Tcs[sortidx]
Tcs_err = Tcs_err[sortidx]
p1 = plot(grid=false, ylims=(0, 1))
plot!(p1, Js, Tcs, color="red", label=nothing, ribbon=Tcs_err)
scatter!(p1, Js, Tcs, color="red", label="LGE Tc")
xlabel!(p1, "J")
# title!(p1, "Tc for V0=$V0, V1=$V1, θ=$(θ_to_π(θ))\n on $size_str $pot lattice")

# Superfluid stiffness 
df_BdG = df_BdG[df_BdG.T.==0, :] # stiffness at 0T only
gdf_BdG = groupby(df_BdG, [:J])
BdG_mean = DataFrame(J=[], Tc=[])

for g in gdf_BdG
    for r in eachrow(g)
        Tc = phase_fluctuation_Tc(m; Δ=r.Δ, K=r.K, Π=r.Π)
        dfi = DataFrame(J=[r.J], Tc=[Tc])
        append!(BdG_mean, dfi)
    end
end

dirs = ["Dₛ/π (x̂)", "Dₛ/π (ŷ)", "Dₛ/π (ẑ)"]
cmap = ["green", "orange", "blue"]
Js, Tcs = BdG_mean.J, BdG_mean.Tc
sortidx = sortperm(Js)
Js = Js[sortidx]
Tcs = Tcs[sortidx]
for i in 1:ndims
    toplot = [Tc[i] for Tc in Tcs]
    plot!(p1, Js, toplot, label=nothing, c=cmap[i], secondary=true)
    scatter!(p1, Js, toplot, label=dirs[i], c=cmap[i], secondary=true)
end
plot!(p1, legend=:right)
plot!(p1, legend=false)

if savefigs
    savefig(p1, joinpath(@__DIR__, "figures", "$(L)L_$(V0)V0_$(V1)V1_$(pot)_stiffness.pdf"))
end

