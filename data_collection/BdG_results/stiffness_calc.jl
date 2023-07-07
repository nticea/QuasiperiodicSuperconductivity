## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")

## PARAMETERS ##
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
periodic = true
T = 0
tol = 1e-15

# load in the data 
savepath_BdG = joinpath(@__DIR__, "BdG_results.csv")
df_BdG = load_dataframe(savepath_BdG)
df_BdG[!, :K] = [zeros(4) for _ in 1:nrow(df_BdG)]
df_BdG[!, :Π] = [zeros(4) for _ in 1:nrow(df_BdG)]

# add column to df for K and Π
for (i, row) in enumerate(eachrow(df_BdG))
    @show i
    @assert row.T == 0
    K, Π = superfluid_stiffness_finiteT(row.T, L=row.L, t=t, J=row.J, Q=Q, μ=μ, V0=row.V0, V1=row.V1, tol=tol, θ=row.θ, ϕx=row.ϕx, ϕy=row.ϕy, niter=1, periodic=periodic, Δ_init=row.Δ)
    df_BdG[i, :K] = K
    df_BdG[i, :Π] = Π
    CSV.write(savepath_BdG, df_BdG)
end


