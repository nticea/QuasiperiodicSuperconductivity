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

## MODEL PARAMETERS ##
args = parse.(Float64, ARGS)
@assert length(args) == 11
L, t, Q, μ, θ, ϕx, ϕy, V0, V1, J, periodic = [args[n] for n in 1:length(args)]
L = Int(L)
periodic = Bool(periodic)

Ts = expspace(-5, 1, 20) # temperature 
symmetries = ["d-wave"] # model symmetry 

for T in Ts
    for symmetry in symmetries

        @show J, symmetry, T

        ## SAVING ##  
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
        mkpath(joinpath(@__DIR__, "data"))
        savepath = joinpath(@__DIR__, "data", "$(L)L_$(J)J_$(symmetry)" * timestamp * ".csv")

        ## Tc using LGE ##
        println("Finding χ")
        χ = @time pairfield_susceptibility(T, symmetry, L=L, t=t, J=J,
            Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic)
        df = DataFrame(L=[L], J=[J], Q=[Q], θ=[θ],
            ϕx=[ϕx], ϕy=[ϕy], symmetry=[symmetry],
            T=[T], χ=[χ])
        CSV.write(savepath, df)
        flush(stdout)

    end
end
