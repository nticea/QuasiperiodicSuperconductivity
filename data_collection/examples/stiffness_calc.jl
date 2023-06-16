## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/stiffness.jl")

## PARAMETERS ##
L = 9 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = nothing#π / 7
ϕx = 0#0.375
ϕy = 0#0.29
V0 = -1
V1 = 0
J = 2
periodic = true
niter = 500
tol = 1e-6
pairing_symmetry = "d-wave"
noise = 0#1e-3 # looks like adding noise is not helping anything 

T = 0.07
superfluid_stiffness(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, noise=noise)