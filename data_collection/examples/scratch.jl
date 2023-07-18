
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")

## PARAMETERS ##
L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
periodic = true
ϕx = 0
ϕy = 0
V0 = -2.3
V1 = 0
J = 1
LGE_tol = 1e-2
BdG_tol = 1e-12
niter = 500

## Tc using LGE ##
println("Finding Tc using LGE")
Tc, λ, Δ_LGE = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, tol=LGE_tol, npts=5)

# Get the corresponding BdG spatial profile 
if isfinite(Tc)
    println("Finding BdG spatial profile at Tc")
    Δ_BdG, hist = compute_Δ_dwave(Tc; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=BdG_tol, Δ_init=Δ_LGE)

    ## Superfluid stiffness calculation ##
    T = 0 # everything is at 0 temperature

    # Get the initial LGE guess 
    println("Finding LGE sol'n at T=0")
    λ, Δ_LGE = @time pairfield_correlation(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)

    # Get the BdG parameters 
    println("Finding BdG coefficients at T=0")
    Δ_BdG, hist = @time compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=BdG_tol, Δ_init=Δ_LGE)

    # Superfluid stiffness
    K, Π = @time superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=BdG_tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=1, periodic=periodic, Δ_init=Δ_BdG)
end
