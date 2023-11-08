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

L = 3 # the full system is L × L 
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
df_LGE_full = load_LGE(dirname)
df_BdG_full = load_BdG(dirname)

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
p1 = plot(grid=false)
plot!(p1, Js, Tcs, color="red", label=nothing, ribbon=Tcs_err)
scatter!(p1, Js, Tcs, color="red", label="LGE Tc")
xlabel!(p1, "J")
title!(p1, "Tc for V0=$V0, V1=$V1, θ=$(θ_to_π(θ))\n on $size_str lattice")

# Superfluid stiffness 
df_BdG = df_BdG[df_BdG.T.==0, :] # stiffness at 0T only
gdf_BdG = groupby(df_BdG, [:J])
BdG_mean = DataFrame(J=[], K_mean=[], Π_mean=[], D_mean=[], K_sem=[], Π_sem=[], D_sem=[])
# combine(gdf_BdG, [:K => mean, :Π => mean,
#     :K => sem, :Π => sem])

function sem_dims(arr; dims)
    @assert Base.ndims(arr) == 2
    res = []
    if dims == 1
        for a in eachcol(arr)
            push!(res, sem(a))
        end
    elseif dims == 2
        for a in eachrow(arr)
            push!(res, sem(a))
        end
    else
        println("idk man :(")
        return
    end
    return res
end

for g in gdf_BdG
    J = g.J[1]
    Ks, Πs = hcat(g.K...), hcat(g.Π...)
    Ks_avg, Πs_avg = mean(Ks, dims=2), mean(Πs, dims=2)
    Ds_avg = -Ks_avg + Πs_avg
    Ks_sem, Πs_sem = sem_dims(Ks, dims=2), sem_dims(Πs, dims=2)
    Ds_sem = -Ks_sem + Πs_sem
    dfi = DataFrame(J=[J], K_mean=[Ks_avg], Π_mean=[Πs_avg], D_mean=[Ds_avg], K_sem=[Ks_sem], Π_sem=[Πs_sem], D_sem=[Ds_sem])
    append!(BdG_mean, dfi)
end

dirs = ["Dₛ/π (x̂)", "Dₛ/π (ŷ)", "Dₛ/π (ẑ)"]
cmap = ["green", "orange", "blue"]
Js, Ds_avg, Ds_sem = BdG_mean.J, BdG_mean.D_mean, BdG_mean.D_sem
sortidx = sortperm(Js)
Js = Js[sortidx]
Ds_avg = Ds_avg[sortidx]
Ds_sem = Ds_sem[sortidx]
for i in 1:ndims
    toplot = [D[i] for D in Ds_avg]
    toplot_ribbon = [D[i] for D in Ds_sem]
    plot!(p1, Js, toplot, label=nothing, c=cmap[i], secondary=true, ribbon=toplot_ribbon)
    scatter!(p1, Js, toplot, label=dirs[i], c=cmap[i], secondary=true)
end
plot!(p1, legend=:right)

if savefigs
    savefig(p1, joinpath(@__DIR__, "figures", "$(L)L_$(V0)V0_$(V1)V1_stiffness_averaged.pdf"))
end

