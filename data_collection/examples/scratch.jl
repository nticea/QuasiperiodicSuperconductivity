
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
L = 11 # the full system is L × L 
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

