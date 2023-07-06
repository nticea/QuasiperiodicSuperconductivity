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
L = 7 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
V0 = 1
V1 = -1.5
periodic = true
niter = 500
tol = 1e-15
T = 0

# saving information 
savepath_BdG = joinpath(@__DIR__, "BdG_results.csv")
savepath_LGE = joinpath(@__DIR__, "LGE_results.csv")

# load in the dataframe, if it exists. If not, make a new one
df_BdG = load_dataframe(savepath_BdG)
df_LGE = load_dataframe(savepath_LGE)

Js = collect(0:0.1:2)
ϕxs = LinRange(0, 2 * π, 5)
ϕys = LinRange(0, 2 * π, 5)

## RUNNING THE CODE ## 
for J in Js # iterate through all J values
    println("Running J=$(J)")
    for ϕx in ϕxs
        for ϕy in ϕys
            if !already_calculated(df_BdG; L=L, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, T=T)
                λ, Δ_LGE = pairfield_correlation(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)
                update_results!(df_LGE; L=L, T=T, λ=λ, Δ=Δ_LGE, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
                CSV.write(savepath_LGE, df_LGE)
                flush(stdout)

                Δ_BdG, hist = compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, Δ_init=Δ_LGE)
                update_results!(df_BdG; L=L, T=T, λ=λ, Δ=Δ_BdG, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
                CSV.write(savepath_BdG, df_BdG)
                flush(stdout)
            end
        end
    end
end


