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
J = 2
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

α = 1
Ts = [0.16, 0.18, 0.2, 0.22]
LGE_sweep(Ts, α=α)

α = 0.8
Ts = [0.2, 0.22, 0.24, 0.28]
LGE_sweep(Ts, α=α)

α = 0.6
Ts = [0.26, 0.3, 0.32, 0.38]
LGE_sweep(Ts, α=α)

α = 0.4
Ts = [0.46, 0.5, 0.55, 0.6, 0.65]
LGE_sweep(Ts, α=α)

α = 0.2
Ts = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35]
LGE_sweep(Ts, α=α)

α = 1.2
Ts = [0.14, 0.16, 0.2, 0.24]
LGE_sweep(Ts, α=α)

α = 1.4
Ts = [0.12, 0.14, 0.16, 0.2]
LGE_sweep(Ts, α=α)

α = 1.6
Ts = [0.16, 0.18, 0.22, 0.24]
LGE_sweep(Ts, α=α)

α = 1.8
Ts = [0.12, 0.16, 0.18, 0.2]
LGE_sweep(Ts, α=α)

α = 2
Ts = [0.12, 0.16, 0.18, 0.2]
LGE_sweep(Ts, α=α)

