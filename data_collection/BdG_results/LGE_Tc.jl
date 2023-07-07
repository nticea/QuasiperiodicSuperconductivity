## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/meanfield.jl")

## PARAMETERS ##
L = 11 # the full system is L × L 
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
savepath_LGE = joinpath(@__DIR__, "LGE_Tc.csv")

# load in the dataframe, if it exists. If not, make a new one
df_LGE = load_dataframe(savepath_LGE)

Js = collect(0:0.1:2)
ϕxs = LinRange(0, π, 5)
ϕys = LinRange(0, π, 5)

## RUNNING THE CODE ## 
for J in Js # iterate through all J values
    println("Running J=$(J)")
    for ϕx in ϕxs
        for ϕy in ϕys
            if !already_calculated(df_LGE; L=L, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
                Tc, λ, Δ_LGE = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)
                update_results!(df_LGE; L=L, T=Tc, λ=λ, Δ=Δ_LGE, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
                CSV.write(savepath_LGE, df_LGE)
                flush(stdout)
            end
        end
    end
end

