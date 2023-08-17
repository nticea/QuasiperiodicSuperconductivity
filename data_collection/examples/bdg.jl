## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

## PARAMETERS ## 
L = 7
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

# simulation parameters
niter = 500
tol = 1e-12

# initialize model 
m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)

# find Δ
Tc, λ, Δ_LGE = LGE_find_Tc(m)
Δ_BdG, hist = compute_Δ(m, T=Tc, niter=niter, tol=tol, Δ_init=Δ_LGE)

# plot things 
plot_spatial_profile(m, Δ=Δ_BdG)
plot_spatial_profile(m, Δ=Δ_LGE)
plot_in_configuration_space(m, Δ=Δ_LGE)
plot_in_configuration_space(m, Δ=Δ_BdG)


