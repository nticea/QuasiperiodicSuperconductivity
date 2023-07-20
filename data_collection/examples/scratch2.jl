
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
L = 17 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
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
T = 0.18 # everything is at 0 temperature
J = 1.1

Mblocks = @time return_M(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)

# Δ = spatial_profile(Δ_LGE, L=L)
# Δ1 = Δ[1, :, :]
# Δ2 = Δ[2, :, :]
# Δ3 = Δ[3, :, :]
# Δ4 = Δ[4, :, :]

# Δ̃1 = circshift(Δ1, (-1, 0)) # mapping between +x̂ and -x̂
# Δ̃2 = circshift(Δ2, (0, 1)) # mapping between +ŷ and -ŷ

# @show maximum(Δ̃1 - Δ3)
# @show maximum(Δ̃2 - Δ4)