## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using CSV
using DataFrames
using StatsPlots

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
V0 = 2#-3
V1 = -1.5#0
ndims = 3
periodic = true
disorder = false
savefigs = false
slice = 1

if disorder
    dirname = "data_$(ndims)D_disorder"
else
    dirname = "data_$(ndims)D_QP"
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
# df_LGE_full = load_LGE(dirname)
# df_BdG_full = load_BdG(dirname)

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
p1 = plot(ylims=(-0.1, 1))
plot!(p1, Js, Tcs, color="red", label=nothing, ribbon=Tcs_err)
scatter!(p1, Js, Tcs, color="red", label="LGE Tc")
xlabel!(p1, "J")
title!(p1, "Tc for V0=$V0, V1=$V1, θ=$(θ_to_π(θ))\n on $size_str lattice")

@assert 1 == 0

# Superfluid stiffness 
df_BdG = df_BdG[df_BdG.T.==0, :] # stiffness at 0T only
gdf_BdG = groupby(df_BdG, [:J])
BdG_mean = DataFrame(["J", "K_mean", "Π_mean", "K_sem", "Π_sem"])
# combine(gdf_BdG, [:K => mean, :Π => mean,
#     :K => sem, :Π => sem])
for 

Js = sort(unique(BdG_mean.J))
Ds_avg = zeros(length(Js), ndims)
Ks_avg = zeros(length(Js), ndims)
Πs_avg = zeros(length(Js), ndims)
for (j, J) in enumerate(Js)
    dfsub = BdG_mean[(df_BdG_0.J.==J), :]
    Ks, Πs = hcat(dfsub.K...), hcat(dfsub.Π...)
    Ks, Πs = mean(Ks, dims=2), mean(Πs, dims=2)
    Ds_avg[j, :] = -Ks + Πs
    Ks_avg[j, :] = Ks
    Πs_avg[j, :] = Πs
end

dirs = ["Dₛ/π (x̂)", "Dₛ/π (ŷ)", "Dₛ/π (ẑ)"]
cmap = ["green", "orange", "blue"]
for i in 1:ndims
    plot!(p1, Js, Ds_avg[:, i], label=nothing, c=cmap[i], secondary=true)
    scatter!(p1, Js, Ds_avg[:, i], label=dirs[i], c=cmap[i], secondary=true)
end
plot!(p1, legend=:right)

if savefigs
    savefig(p1, joinpath(@__DIR__, "figures", "$(L)L_$(V0)V0_$(V1)V1_stiffness_averaged.pdf"))
end

