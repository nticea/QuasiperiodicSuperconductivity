## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/stiffness.jl")

## PARAMETERS ## 
L = 17
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

# Simulation parameters
niter = 500
tol = 1e-12

# Find LGE soln at T=0 
λ, ΔLGE = pairfield_correlation(m, T=0)

# Calculate superfluid stiffness
K, Π, Δ = superfluid_stiffness_finiteT(m, T=0, niter=niter, tol=tol, Δ_init=ΔLGE)


