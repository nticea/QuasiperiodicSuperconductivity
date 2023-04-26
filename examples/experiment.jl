## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using LinearAlgebra, Arpack, Plots
include("../src/model.jl")
include("../src/meanfield.jl")

## PARAMETERS ##
L = 10 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = nothing
J = 2
T = 0
V0 = 1
V1 = 0
periodic = true

λ = λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
