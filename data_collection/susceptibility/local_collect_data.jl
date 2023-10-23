## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames
using Dates

include("../../src/meanfield.jl")
include("utilities.jl")

L = 7
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
periodic = true # periodic 
disorder = false
ndims = 3
nrep = 40
Js = [1, 4]
Ts = expspace(-3, 1, 30) # temperature 

Λ = nothing
if disorder
    dirname = "data_random"
else
    dirname = "data"
end
if periodic
    dirname = dirname * "_PBC"
else
    dirname = dirname * "_OBC"
end

# datapath
datapath = joinpath(@__DIR__, dirname)
mkpath(datapath)

for _ in 1:nrep
    ϕx = 2π * rand()
    ϕy = 2π * rand()
    ϕz = 2π * rand()

    for J in Js
        for T in Ts
            m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)

            ## CALCULATION ## 
            χ, dχdlogT = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true, Λ=Λ)

            ## SAVING ##  
            timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
            savepath = joinpath(datapath, "$(L)L_$(J)J" * timestamp * ".csv")
            df = DataFrame(L=[L], J=[J], Q=[Q], θ=[θ],
                ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ndims=[ndims],
                T=[T], χ=[χ], dχdlogT=[dχdlogT], Λ=[6])
            CSV.write(savepath, df)
        end
    end
end

