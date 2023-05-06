## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using LinearAlgebra, Arpack, Plots
using Profile
include("../src/model.jl")
include("../src/meanfield.jl")
# include("../src/BdG_dwave.jl")
# include("../src/BdG.jl")

## PARAMETERS ##
L = 30 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = nothing
J = 3
T = 0.24
V0 = 1
V1 = -1.5
periodic = true

@time M = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="s-wave")
@time M = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")

println("Performing eigendecomposition")

@time decomp, _ = partialschur(Hermitian(M), nev=1, tol=1e-6, which=LM())
maxev = decomp.Q
λ = decomp.R[1]
@show λ