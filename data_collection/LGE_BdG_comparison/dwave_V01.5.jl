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
include("../../src/BdG_dwave.jl")

## PARAMETERS ##
L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0.375
ϕy = 0.29
V0 = 1.5
V1 = -1
niter = 500
tol = 1e-12
periodic = true
pairing_symmetry = "d-wave"
num_T = 20

# saving information 
savepath_LGE = joinpath(@__DIR__, "data", "$(L)L$(V0)V0$(V1)V1_LGE_results.csv")
savepath_BdG = joinpath(@__DIR__, "data", "$(L)L$(V0)V0$(V1)V1_BdG_results.csv")

# load in the dataframe, if it exists. If not, make a new one
df_LGE = load_dataframe(savepath_LGE)
df_BdG = load_dataframe(savepath_BdG)

Js = collect(0:0.1:2)
for J in Js
    # find Tc and save LGE result 
    if !already_calculated(df_BdG; L=L, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
        Tc, λ, Δ_LGE = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)
        update_results!(df_LGE; L=L, T=Tc, λ=λ, Δ=Δ_LGE, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
        CSV.write(savepath_LGE, df_LGE)

        if abs(λ - 1) < 1e-3
            # run some BdG calculations in the range of this Tc 
            Ts = LinRange(Tc - 0.1, Tc + 0.1, num_T)
            for T in Ts
                Δ_BdG, hist = compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, Δ_init=Δ_LGE)
                update_results!(df_BdG; L=L, T=T, λ=0, Δ=Δ_BdG, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1)
                CSV.write(savepath_BdG, df_BdG)
                flush(stdout)
            end
        end
    end
end

