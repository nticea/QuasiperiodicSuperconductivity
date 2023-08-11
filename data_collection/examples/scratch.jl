
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames
using CurveFit

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")

L = 17
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0
ϕy = 0
periodic = true
J = 2
T = 0

stamp = "diagonalized_$(L)L_$(J)J_$(round(θ, digits=3))theta_$(round(Q,digits=3))Q.h5"
scratchpath = joinpath("/scratch/users/nticea", "QuasiperiodicSuperconductivity", "diagonalized_hamiltonians", stamp)

H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
E, U = diagonalize_hamiltonian(H0)
DH = DiagonalizedHamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic, E=E, U=U)
save_structs(DH, stamp)

uniform_susceptibility(T, L=L, t=t, J=J,
    Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic)

# times = []
# bytes = []
# for L in Ls
#     @show L
#     H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
#     # Diagonalize this Hamiltonian
#     res = @timed diagonalize_hamiltonian(H0)

#     push!(times, res.time)
#     push!(bytes, res.bytes)
# end
