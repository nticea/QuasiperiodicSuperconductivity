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
J = 3
T = 1e-1
V0 = 0.2
V1 = 0
periodic = true

M, λ = λmax(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, symmetry="d-wave")
heatmap(reverse(M, dims=2), c=:bwr, clims=(-maximum(abs.(M)), maximum(abs.(M))))
@show λ