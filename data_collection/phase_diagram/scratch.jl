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
include("utilities.jl")

## MODEL PARAMETERS ##
L = 17 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0
ϕy = 0
periodic = true
V1 = -2
V0 = -0.5
J = 1.5

## SIMULATION PARAMETERS ## 
niter = 500
BdG_tol = 1e-15
LGE_tol = 1e-2

# check if this point has already been computed 
if !(already_computed(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, J=J, periodic=periodic))
    ## Tc using LGE ##
    println("Finding Tc using LGE")
    Tc, λ, Δ_LGE = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, tol=LGE_tol, npts=5)
end
