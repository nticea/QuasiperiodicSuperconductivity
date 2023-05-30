## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/model.jl")
include("../../src/meanfield.jl")

## PARAMETERS ##
L = 55 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
J = 3
V0 = 1
V1 = -1.5
periodic = true
pairing_symmetry = "d-wave"

# temperature range 
Ts = [0.24, 0.23, 0.22]

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_results.csv")

# load in the dataframe, if it exists. If not, make a new one
df = load_dataframe(savepath)

## RUNNING THE CODE ## 
println("Running J=$(J)")
iter = ProgressBar(1:length(Ts))
for i in iter # iterate through all temperatures
    T = Ts[i]
    # check whether this particular (J,T,V0) combo has been already computed 
    if !already_calculated(df; L=L, J=J, θ=θ, V0=V0, V1=V1, T=T)

        # calculate M at a given (J,T,V0, V1)
        @time λ, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")

        update_results!(df; L=L, λ=λ, J=J, θ=θ, V0=V0, V1=V1, T=T, Δ=Δ)
        CSV.write(savepath, df)
        flush(stdout)
    end
end
