
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates
include("../../src/meanfield.jl")

## PARAMETERS ## 
L = 7
J = 1
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
ϕx = 0
ϕy = 0
periodic = true
ndims = 3

g = Graphs.SimpleGraphs.grid((L, L, L), periodic=periodic)
H = Graphs.LinAlg.adjacency_matrix(g)

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
