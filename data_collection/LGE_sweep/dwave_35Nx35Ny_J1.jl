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
L = 35 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
J = 1
V1 = -1
periodic = true
pairing_symmetry = "d-wave"

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_$(J)J_dwave_results.csv")

# load in the dataframe, if it exists. If not, make a new one
df = load_dataframe(savepath)

function LGE_sweep(Ts; α::Real)
    ## RUNNING THE CODE ## 
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))

    # Compute V1, holding the ratio α≡V0/V1 fixed 
    V0 = V1 * α

    for i in iter # iterate through all temperatures
        T = Ts[i]
        # check whether this particular (J,T,V0) combo has been already computed 
        if !already_calculated(df; L=L, J=J, θ=θ, V0=V0, V1=V1, T=T)

            # calculate M at a given (J,T,V0, V1)
            @time λ, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)

            update_results!(df; L=L, λ=λ, J=J, θ=θ, V0=V0, V1=V1, T=T, Δ=Δ)
            CSV.write(savepath, df)
            flush(stdout)
        end
    end
end

α = 0.3
Ts = [0.1, 0.12, 0.15, 0.17]
LGE_sweep(Ts, α=α)

α = 0.425
LGE_sweep(Ts, α=α)

α = 0.55
LGE_sweep(Ts, α=α)

α = 0.675
LGE_sweep(Ts, α=α)

α = 0.8
LGE_sweep(Ts, α=α)

α = 0.925
LGE_sweep(Ts, α=α)

α = 1.05
Ts = [0.12, 0.15, 0.17, 0.19]
LGE_sweep(Ts, α=α)

α = 1.175
Ts = [0.12, 0.15, 0.17, 0.19]
LGE_sweep(Ts, α=α)

α = 1.3
Ts = [0.15, 0.17, 0.19, 0.22]
LGE_sweep(Ts, α=α)

