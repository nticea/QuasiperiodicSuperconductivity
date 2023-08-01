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

Ts = expspace(-4, 1, 50) # temperature 

for T in Ts
    @show L, J, T

    ## CALCULATION ## 
    println("Finding χ")
    χ = @time uniform_susceptibility(T, L=L, t=t, J=J,
        Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic)

    ## SAVING ##  
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
    mkpath(joinpath(@__DIR__, "data"))
    savepath = joinpath(@__DIR__, "data", "$(L)L_$(J)J" * timestamp * ".csv")
    df = DataFrame(L=[L], J=[J], Q=[Q], θ=[θ],
        ϕx=[ϕx], ϕy=[ϕy],
        T=[T], χ=[χ])
    CSV.write(savepath, df)
    flush(stdout)
end
