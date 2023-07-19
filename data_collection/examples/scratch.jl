
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
V0 = 1
V1 = -1.5
LGE_tol = 1e-2
BdG_tol = 1e-12
niter = 500

## Tc using LGE ##

## Superfluid stiffness calculation ##
T = 0 # everything is at 0 temperature

Ks, Πs = [], []
Js = [0, 1, 2, 3, 4]
for J in Js
    # Get the initial LGE guess 
    println("Finding LGE sol'n at T=0")
    λ, Δ_LGE = @time pairfield_correlation(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)

    # Δ_BdG, hist = compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=BdG_tol, Δ_init=Δ_LGE)
    # Get the BdG parameters 
    println("Finding BdG coefficients at T=0")

    # Superfluid stiffness
    K, Π, Δ = @time superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=BdG_tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, Δ_init=Δ_LGE)
    push!(Ks, K[1])
    push!(Πs, Π[1])
end
# Δ0 = zeros(size(Δ_LGE)...)
# K, Π, Δ = @time superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=BdG_tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=1, periodic=periodic, Δ_init=Δ0)

plot(Js, Ks, label="K")
plot!(Js, Πs, label="Π")

plot(Js, -Ks + Πs, label="Dₛ/π")