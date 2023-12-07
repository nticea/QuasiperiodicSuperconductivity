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
    dirname = "FINAL_$(ndims)D_data_random"
    println("We are adding disorder!")
else
    dirname = "FINAL_$(ndims)D_data"
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

# all the other things we have computed 
Js = collect(0:0.05:1.25)
Js = shuffle(Js)

for J in Js
    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
    @show L, J, T, ndims

    ## CALCULATION ## 
    println("Finding χ")
    χ0_dwave, dχdlogT_dwave, χ0_pwave, dχdlogT_pwave = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true, Λ=Λ)

    ## SAVING ##  
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
    savepath = joinpath(datapath, "$(L)L_$(J)J" * timestamp * ".csv")
    df = DataFrame(L=[L], J=[J], Q=[Q], θ=[θ], μ=[μ],
        ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ndims=[ndims], T=[T], Λ=[Λ],
        χ_dwave=[χ0_dwave], dχdlogT_dwave=[dχdlogT_dwave],
        χ_pwave=[χ0_pwave], dχdlogT_pwave=[dχdlogT_pwave])
    CSV.write(savepath, df)
    flush(stdout)
end

