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
L = 17 # the full system is L × L 
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

# saving information 
savepath_LGE = joinpath(@__DIR__, "data", "$(L)L$(V0)V0$(V1)V1_LGE_results.csv")
savepath_BdG = joinpath(@__DIR__, "data", "$(L)L$(V0)V0$(V1)V1_BdG_results.csv")

# load in the dataframe, if it exists. If not, make a new one
df_LGE = load_dataframe(savepath_LGE)
df_BdG = load_dataframe(savepath_BdG)

function LGE_sweep(Ts; J::Real)
    ## RUNNING THE CODE ## 
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))

    for i in iter # iterate through all temperatures
        T = Ts[i]
        # check whether this particular (J,T,V0) combo has been already computed 
        if !already_calculated(df_LGE; L=L, J=J, θ=θ, V0=V0, V1=V1, T=T, ϕx=ϕx, ϕy=ϕy)

            # calculate the LGE result 
            @time λ, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)

            # save the LGE results 
            update_results!(df_LGE; L=L, λ=λ, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, T=T, Δ=Δ_LGE)
            CSV.write(savepath_LGE, df_LGE)
            flush(stdout)

            if !already_calculated(df_BdG; L=L, J=J, θ=θ, V0=V0, V1=V1, T=T, ϕx=ϕx, ϕy=ϕy)

                # calculate the BdG result 
                @time Δ_BdG, BdG_hist = compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, Δ_init=Δ_LGE)

                # save the BdG results 
                update_results!(df_BdG; L=L, λ=λ, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, T=T, Δ=Δ_BdG)
                CSV.write(savepath_BdG, df_BdG)
                flush(stdout)
            end
        end
    end
end

J = 0
Ts = LinRange(0.12, 0.16, 10)
LGE_sweep(Ts, J=J)

J = 0.1
Ts = LinRange(0.12, 0.16, 10)
LGE_sweep(Ts, J=J)

J = 0.2
Ts = LinRange(0.11, 0.15, 10)
LGE_sweep(Ts, J=J)

J = 0.3
Ts = LinRange(0.10, 0.14, 10)
LGE_sweep(Ts, J=J)

J = 0.4
Ts = LinRange(0.08, 0.12, 10)
LGE_sweep(Ts, J=J)

J = 0.5
Ts = LinRange(0.07, 0.11, 10)
LGE_sweep(Ts, J=J)

J = 0.6
Ts = LinRange(0.07, 0.11, 10)
LGE_sweep(Ts, J=J)

J = 0.7
Ts = LinRange(0.08, 0.12, 10)
LGE_sweep(Ts, J=J)

J = 0.8
Ts = LinRange(0.09, 0.13, 10)
LGE_sweep(Ts, J=J)

J = 0.9
Ts = LinRange(0.1, 0.14, 10)
LGE_sweep(Ts, J=J)

J = 1
Ts = LinRange(0.11, 0.15, 10)
LGE_sweep(Ts, J=J)

J = 1.1
Ts = LinRange(0.11, 0.15, 10)
LGE_sweep(Ts, J=J)

J = 1.2
Ts = LinRange(0.12, 0.16, 10)
LGE_sweep(Ts, J=J)

J = 1.3
Ts = LinRange(0.125, 0.165, 10)
LGE_sweep(Ts, J=J)

J = 1.4
Ts = LinRange(0.13, 0.17, 10)
LGE_sweep(Ts, J=J)

J = 1.5
Ts = LinRange(0.13, 0.17, 10)
LGE_sweep(Ts, J=J)

J = 1.6
Ts = LinRange(0.135, 0.175, 10)
LGE_sweep(Ts, J=J)

J = 1.7
Ts = LinRange(0.14, 0.18, 10)
LGE_sweep(Ts, J=J)

J = 1.8
Ts = LinRange(0.145, 0.185, 10)
LGE_sweep(Ts, J=J)

J = 1.9
Ts = LinRange(0.15, 0.19, 10)
LGE_sweep(Ts, J=J)

J = 2
Ts = LinRange(0.1525, 0.1925, 10)
LGE_sweep(Ts, J=J)
