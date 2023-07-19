
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
L = 7 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
periodic = true
ϕx = 0
ϕy = 0
V0 = -2.3
V1 = 0
J = 0
LGE_tol = 1e-2
BdG_tol = 1e-12
niter = 500

## Tc using LGE ##

## Superfluid stiffness calculation ##
T = 0 # everything is at 0 temperature

# Get the initial LGE guess 
println("Finding LGE sol'n at T=0")
λ, Δ_LGE = @time pairfield_correlation(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)

# Get the BdG parameters 
println("Finding BdG coefficients at T=0")

# Superfluid stiffness
K, Π, Δ = @time superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=BdG_tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, Δ_init=Δ_LGE)

