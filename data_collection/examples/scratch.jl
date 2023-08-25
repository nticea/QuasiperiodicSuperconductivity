## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

## PARAMETERS ## 
L = 23
J = 0.8
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 1.5
V1 = -1
ϕx = 0
ϕy = 0
ϕz = 0
periodic = true
ndims = 3

# initialize model 
m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
E, U = @time diagonalize_hamiltonian(m)
plot(E)

histogram(E, bins=100, xlabel="E", ylabel="Frequency", title="Histogram")
