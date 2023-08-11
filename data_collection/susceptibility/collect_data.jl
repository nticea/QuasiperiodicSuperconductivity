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

## MODEL PARAMETERS ##
args = parse.(Float64, ARGS)
@assert length(args) == 11
L, t, Q, μ, θ, ϕx, ϕy, V0, V1, J, periodic = [args[n] for n in 1:length(args)]
L = Int(L)
periodic = Bool(periodic)

Ts = expspace(-3, 1, 30) # temperature 
Λ = 0.3

# all the other things we have computed 
dfs = load_dfs()

# diagonalize and save the Hamiltonian 
stamp = "diagonalized_$(L)L_$(J)J_$(round(θ, digits=3))theta_$(round(Q,digits=3))Q.h5"
scratchbase = joinpath("/scratch/users/nticea", "QuasiperiodicSuperconductivity", "diagonalized_hamiltonians")
mkdir(scratchbase)
scratchpath = joinpath(scratchbase, stamp)

H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
E, U = diagonalize_hamiltonian(H0)
DH = DiagonalizedHamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic, E=E, U=U)
save_structs(DH, scratchpath)

for T in Ts
    if !already_computed(dfs, T=T, L=L, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, J=J, Λ=Λ)
        @show L, J, T

        ## CALCULATION ## 
        println("Finding χ")
        χ = @time uniform_susceptibility(T, L=L, t=t, J=J,
            Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic,
            checkpointpath=scratchpath)

        ## SAVING ##  
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
        mkpath(joinpath(@__DIR__, "data"))
        savepath = joinpath(@__DIR__, "data", "$(L)L_$(J)J" * timestamp * ".csv")
        df = DataFrame(L=[L], J=[J], Q=[Q], θ=[θ],
            ϕx=[ϕx], ϕy=[ϕy],
            T=[T], χ=[χ], Λ=[Λ])
        CSV.write(savepath, df)
        flush(stdout)
    end
end
