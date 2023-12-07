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

L = 11
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0.75
θ = π / 7
V0 = 0
V1 = 0
ndims = 3
disorder = false
periodic = true

Js = collect(0:0.05:1.25)
Ts = expspace(-3, -1, 15) # temperature 

nrep = 40

# datapath
disorder_datapath = joinpath(@__DIR__, "data_random_$(ndims)D")
QP_datapath = joinpath(@__DIR__, "data_QP_$(ndims)D")
mkpath(disorder_datapath)
mkpath(QP_datapath)

Λ = 0.3
for n in 1:nrep
    ϕx = 2π * rand()
    ϕy = 2π * rand()
    ϕz = 2π * rand()

    for J in Js
        for T in Ts
            @show n
            ## QUASIPERIODIC 
            m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims, disorder=false)
            χ0_dwave, dχdlogT_dwave, χ0_pwave, dχdlogT_pwave = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true, Λ=Λ)
            ## SAVING ##  
            timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
            savepath = joinpath(QP_datapath, "$(L)L_$(J)J" * timestamp * ".csv")
            df = DataFrame(L=[L], J=[J], Q=[Q], θ=[θ], μ=[μ],
                ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ndims=[ndims], T=[T], Λ=[Λ],
                χ_dwave=[χ0_dwave], dχdlogT_dwave=[dχdlogT_dwave],
                χ_pwave=[χ0_pwave], dχdlogT_pwave=[dχdlogT_pwave])
            CSV.write(savepath, df)

            ## DISORDER 
            m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims, disorder=true)
            χ0_dwave, dχdlogT_dwave, χ0_pwave, dχdlogT_pwave = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true, Λ=Λ)
            timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
            savepath = joinpath(disorder_datapath, "$(L)L_$(J)J" * timestamp * ".csv")
            df = DataFrame(L=[L], J=[J], Q=[Q], θ=[θ], μ=[μ],
                ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ndims=[ndims], T=[T], Λ=[Λ],
                χ_dwave=[χ0_dwave], dχdlogT_dwave=[dχdlogT_dwave],
                χ_pwave=[χ0_pwave], dχdlogT_pwave=[dχdlogT_pwave])
            CSV.write(savepath, df)
        end
    end
end


