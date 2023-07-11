## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames
using Dates

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")

## DELETE THIS ##
L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0 # fixed 
ϕy = 0 # fixed 
periodic = 1

# Parameters to iterate through 
V0 = 1 # iterate through this 
J = 2 # iterate through this 
V1 = -1.5 # look at just a few values: V1=0, -0.5, -1, -1.5

## MODEL PARAMETERS ##
args = parse.(Float64, ARGS)
@assert length(args) == 11
L, t, Q, μ, θ, ϕx, ϕy, V0, V1, J, periodic = [args[n] for n in 1:length(args)]
L = Int(L)
periodic = Bool(periodic)

## SIMULATION PARAMETERS ## 
niter = 500
BdG_tol = 1e-15
LGE_tol = 1e-2

## SAVING ##  
timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
savepath_LGE = joinpath(@__DIR__, "data", "$(L)L_ΦQ_LGE_" * timestamp * ".csv")
df_LGE = load_dataframe(savepath_LGE)

## Tc using LGE ##
println("Finding Tc using LGE")
Tc, λ, Δ_LGE = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, tol=LGE_tol, npts=5)
update_results!(df_LGE; L=L, T=Tc, λ=λ, Δ=Δ_LGE, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
CSV.write(savepath_LGE, df_LGE)
flush(stdout)