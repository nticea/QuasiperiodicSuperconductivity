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
## MODEL PARAMETERS ##
args = parse.(Float64, ARGS)
@assert length(args) == 16
L, t, Q, μ, θ, ϕx, ϕy, ϕz, V0, V1, J, periodic, ndims, disorder, λ, E₀ = [args[n] for n in 1:length(args)]
L = Int(L)
ndims = Int(ndims)
periodic = Bool(periodic)
disorder = Bool(disorder)

datapath = joinpath(@__DIR__, "data")
mkpath(datapath)

Js = LinRange(0, 5, 30)
for J in Js
    @show L, J

    # # for saving hamiltonian
    # stamp = "diagonalized_$(ndims)D_$(L)L_$(J)J_$(round(θ, digits=3))theta_$(round(Q,digits=3))Q_$(round(ϕx,digits=3))ϕx_$(round(ϕy,digits=3))ϕy_$(round(ϕz,digits=3))ϕz_$(Int(periodic))periodic.h5"
    # scratchbase = joinpath("/scratch/users/nticea", "QuasiperiodicSuperconductivity", "diagonalized_hamiltonians")
    # #scratchbase = joinpath(@__DIR__, "diagonalized_hamiltonians")
    # mkpath(scratchbase)
    # scratchpath = joinpath(scratchbase, stamp)

    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
    α₀, ipr_real, ipr_k = compute_scaling_properties(m; λ=λ, E₀=E₀)
    @show α₀, ipr_real, ipr_k
    # save it out 
    save_results(m, λ=λ, E₀=E₀, α₀=α₀, ipr_real=ipr_real, ipr_k=ipr_k)
end