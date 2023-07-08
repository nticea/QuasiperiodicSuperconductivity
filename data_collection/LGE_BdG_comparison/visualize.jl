## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/results.jl")
include("../../src/meanfield.jl")
include("../../src/BdG_dwave.jl")

## PARAMETERS 

L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0.375
ϕy = 0.29
V0 = 1.5
V1 = -1

J = 1.6

## SCRIPT 

df_BdG = load_dataframe("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_BdG_comparison/data/$(L)L$(V0)V0$(V1)V1_BdG_results.csv")
df_LGE = load_dataframe("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_BdG_comparison/data/$(L)L$(V0)V0$(V1)V1_LGE_results.csv")

sub_BdG = df_BdG[(df_BdG.L.==L).&(df_BdG.J.==J).&(df_BdG.V0.==V0).&(df_BdG.V1.==V1).&(df_BdG.θ.==θ).&(df_BdG.ϕx.==ϕx).&(df_BdG.ϕy.==ϕy), :]
sub_LGE = df_LGE[(df_LGE.L.==L).&(df_LGE.J.==J).&(df_LGE.V0.==V0).&(df_LGE.V1.==V1).&(df_LGE.θ.==θ).&(df_LGE.ϕx.==ϕx).&(df_LGE.ϕy.==ϕy), :]

Tc = sub_LGE.T[1]
Ts, Δs = sub_BdG.T, sub_BdG.Δ
norm_Δ = []
for Δ in Δs
    push!(norm_Δ, norm(Δ))
end
p = plot()
plot!(p, Ts, norm_Δ, color="blue", label="BdG", yaxis=:log10)
scatter!(p, Ts, norm_Δ, color="blue", label=nothing, yaxis=:log10)
vline!(p, [Tc], color="black", label="LGE Tc", yaxis=:log10)
title!(p, "BdG Δ(J=$J, V0=$V0, V1=$V1, θ=$(θ_to_π(θ)), \n ϕx=$ϕx,ϕy=$ϕy), for $L × $L lattice")
xlabel!(p, "T")
ylabel!("||Δ||")