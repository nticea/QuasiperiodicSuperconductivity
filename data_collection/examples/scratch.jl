
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/meanfield.jl")

## PARAMETERS ##
L = 17 # the full system is L × L 
t = 1
J = 6.5
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
ϕx = 0
ϕy = 0
V0 = 1
V1 = -1.5
periodic = true

pairfield_correlation(0.2; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)
