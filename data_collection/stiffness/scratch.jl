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

L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = -1
V1 = 0
periodic = true
ϕx = 0
ϕy = 0
J = 1

## SIMULATION PARAMETERS ## 
niter = 500
BdG_tol = 1e-15
LGE_tol = 1e-2

## SAVING ##  
timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
savepath_BdG = joinpath(@__DIR__, "data", "$(L)L_ΦQ_BdG_" * timestamp * ".csv")
savepath_LGE = joinpath(@__DIR__, "data", "$(L)L_ΦQ_LGE_" * timestamp * ".csv")
df_BdG = load_dataframe(savepath_BdG)
df_LGE = load_dataframe(savepath_LGE)

## Tc using LGE ##
println("Finding Tc using LGE")
Tc, λ, Δ_LGE = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, tol=LGE_tol, npts=5)
update_results!(df_LGE; L=L, T=Tc, λ=λ, Δ=Δ_LGE, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
CSV.write(savepath_LGE, df_LGE)
flush(stdout)

# Get the corresponding BdG spatial profile 
if isfinite(Tc)
    ## Superfluid stiffness calculation ##
    T = 0 # everything is at 0 temperature

    # Get the BdG parameters 
    println("Finding BdG coefficients at T=0")
    Δ_BdG, hist = @time compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=BdG_tol)
    flush(stdout)

    # Superfluid stiffness
    K, Π = @time superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=BdG_tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=1, periodic=periodic, Δ_init=Δ_BdG)
    update_results!(df_BdG; L=L, T=T, λ=λ, Δ=Δ_BdG, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, K=K, Π=Π)
    CSV.write(savepath_BdG, df_BdG)
    flush(stdout)
end