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
L, t, Q, μ, θ, _, _, _, V0, V1, J, periodic, ndims, ℓ, E₀ = [args[n] for n in 1:length(args)]
L = Int(L)
ndims = Int(ndims)
periodic = Bool(periodic)

datapath = joinpath(@__DIR__, "data")
mkpath(datapath)

ϕxs = LinRange(0, π + 0.01, 3)
ϕys = LinRange(0, π + 0.1, 3)
ϕzs = LinRange(0, π - 0.01, 3)

ϕBx = 0
ϕBy = 0
ϕBz = 0

dfs = load_dfs()

for ϕx in ϕxs
    for ϕy in ϕys
        for ϕz in ϕzs
            @show L, J

            # for saving hamiltonian
            stamp = "diagonalized_$(ndims)D_$(L)L_$(J)J_$(round(θ, digits=3))theta_$(round(Q,digits=3))Q_$(round(ϕx,digits=3))ϕx_$(round(ϕy,digits=3))ϕy_$(round(ϕz,digits=3))ϕz_$(Int(periodic))periodic.h5"
            scratchbase = joinpath("/scratch/users/nticea", "QuasiperiodicSuperconductivity", "diagonalized_hamiltonians")
            #scratchbase = joinpath(@__DIR__, "diagonalized_hamiltonians")
            mkpath(scratchbase)
            scratchpath = joinpath(scratchbase, stamp)

            m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, ϕBx=ϕBx, ϕBy=ϕBy, ϕBz=ϕBz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)
            if !already_computed(dfs, m; E₀=E₀, ℓ=ℓ)
                α₀, ipr_real, ipr_k = compute_scaling_properties(m; loadpath=scratchpath)
                @show α₀, ipr_real, ipr_k
                # save it out 
                save_results(m, ℓ=ℓ, E₀=E₀, α₀=α₀, ipr_real=ipr_real, ipr_k=ipr_k)
            end
        end
    end
end