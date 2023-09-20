## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using CSV
using DataFrames
using Dates

include("../../src/model.jl")
include("utilities.jl")

## MODEL PARAMETERS ##
args = parse.(Float64, ARGS)
@assert length(args) == 15
L, t, Q, μ, θ, ϕx, ϕy, ϕz, V0, V1, J, periodic, ndims, ℓ, E₀ = [args[n] for n in 1:length(args)]
L = Int(L)
ndims = Int(ndims)
periodic = Bool(periodic)

# for saving hamiltonian
stamp = "diagonalized_$(ndims)D_$(L)L_$(J)J_$(round(θ, digits=3))theta_$(round(Q,digits=3))Q_$(round(ϕx,digits=3))ϕx_$(round(ϕy,digits=3))ϕy_$(round(ϕz,digits=3))ϕz.h5"
scratchbase = joinpath("/scratch/users/nticea", "QuasiperiodicSuperconductivity", "diagonalized_hamiltonians")
# scratchbase = joinpath(@__DIR__, "diagonalized_hamiltonians")
mkpath(scratchbase)
scratchpath = joinpath(scratchbase, stamp)

## DELETE THIS ## 
## END DELETE ##

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)
@time α₀ = multifractal_mean(m; E₀=E₀, ℓ=ℓ, loadpath=scratchpath)

# save it out 
save_results(m, ℓ=ℓ, E₀=E₀, α₀=α₀)