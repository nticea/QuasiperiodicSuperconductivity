## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames
using Dates, Random

include("../../src/meanfield.jl")
include("utilities.jl")

## MODEL PARAMETERS ##
args = parse.(Float64, ARGS)
@assert length(args) == 15
L, t, Q, μ, θ, ϕx, ϕy, ϕz, V0, V1, J, periodic, ndims, disorder, T = [args[n] for n in 1:length(args)]
L = Int(L)
ndims = Int(ndims)
periodic = Bool(periodic)
disorder = Bool(disorder)

Λ = 0.3
if disorder
    dirname = "$(ndims)D_data_random"
    println("We are adding disorder!")
else
    dirname = "$(ndims)D_data"
end
if periodic
    dirname = dirname * "_PBC"
else
    dirname = dirname * "_OBC"
end

# datapath
datapath = joinpath(@__DIR__, dirname)
mkpath(datapath)
# for saving the Hamiltonian
stamp = "diagonalized_$(ndims)D_$(L)L_$(J)J_$(round(θ, digits=3))theta_$(round(Q,digits=3))Q.h5"
scratchbase = joinpath("/scratch/users/nticea", "QuasiperiodicSuperconductivity", "diagonalized_hamiltonians")
# scratchbase = @__DIR__
mkpath(scratchbase)
scratchpath = joinpath(scratchbase, stamp)

# all the other things we have computed 
Js = [0, 0.1, 0.25, 0.3, 0.4, 0.45, 0.5, 0.55]
Js = shuffle(Js)

for J in Js
    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
    @show L, J, T, ndims

    ## CALCULATION ## 
    println("Finding χ")
    χ, dχdlogT = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true, Λ=Λ)

    ## SAVING ##  
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
    savepath = joinpath(datapath, "$(L)L_$(J)J" * timestamp * ".csv")
    df = DataFrame(L=[L], J=[J], Q=[Q], θ=[θ], μ=[μ],
        ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ndims=[ndims],
        T=[T], χ=[χ], dχdlogT=[dχdlogT], Λ=[Λ])
    CSV.write(savepath, df)
    flush(stdout)
end

