## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
include("../src/BdG.jl")

using Plots
using ProgressBars
using CSV
using DataFrames

## PARAMETERS ##
L = 55 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
J = 3
T = 0
V0 = 0.4
pairing_symmetry = "s-wave"
tol = 1e-3

# rund the BdG code 
Δ = compute_Δ(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, θ=θ, tol=tol)