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
Q = 1 / 2 * (√5 - 1)
μ = 1e-8
θ = π / 7
ϕx = 0.375
ϕy = 0.29
V0 = 1.5
V1 = -1
periodic = true
pairing_symmetry = "d-wave"

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_dwave_results.csv")

# load in the dataframe, if it exists. If not, make a new one
df = load_dataframe(savepath)

function LGE_sweep(Ts; J::Real)
    ## RUNNING THE CODE ## 
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))

    for i in iter # iterate through all temperatures
        T = Ts[i]
        # check whether this particular (J,T,V0) combo has been already computed 
        if !already_calculated(df; L=L, J=J, θ=θ, V0=V0, V1=V1, T=T, ϕx=ϕx, ϕy=ϕy)

            # calculate M at a given (J,T,V0, V1)
            @time λ, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)

            update_results!(df; L=L, λ=λ, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, T=T, Δ=Δ)
            CSV.write(savepath, df)
            flush(stdout)
        end
    end
end

# L = 11
# J = 2
# T = 10^(-1.15)
# λ, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")

# @assert 1 == 0

J = 0
Ts = expspace(-0.9, -0.95, 5)
LGE_sweep(Ts, J=J)

J = 0.1
Ts = expspace(-0.8, -0.9, 5)
LGE_sweep(Ts, J=J)

J = 0.2
Ts = expspace(-0.85, -0.95, 5)
LGE_sweep(Ts, J=J)

J = 0.3
Ts = expspace(-0.85, -0.9, 5)
LGE_sweep(Ts, J=J)

J = 0.4
Ts = expspace(-0.85, -0.95, 5)
LGE_sweep(Ts, J=J)

J = 0.5
Ts = expspace(-0.88, -0.92, 5)
LGE_sweep(Ts, J=J)

J = 0.6
Ts = expspace(-0.9, -0.95, 5)
LGE_sweep(Ts, J=J)

J = 0.7
Ts = expspace(-0.95, -1, 5)
LGE_sweep(Ts, J=J)

J = 0.8
Ts = expspace(-0.98, -1.05, 5)
LGE_sweep(Ts, J=J)

J = 0.9
Ts = expspace(-1.05, -1.1, 5)
LGE_sweep(Ts, J=J)

J = 1
Ts = expspace(-1.15, -1.1, 5)
LGE_sweep(Ts, J=J)

J = 1.1
Ts = expspace(-1.17, -1.23, 5)
LGE_sweep(Ts, J=J)

J = 1.2
Ts = expspace(-1.28, -1.32, 5)
LGE_sweep(Ts, J=J)

J = 1.3
Ts = expspace(-1.32, -1.37, 5)
LGE_sweep(Ts, J=J)

J = 1.4
Ts = expspace(-1.32, -1.37, 5)
LGE_sweep(Ts, J=J)

J = 1.5
Ts = expspace(-1.3, -1.35, 5)
LGE_sweep(Ts, J=J)

J = 1.6
Ts = expspace(-1.25, -1.3, 5)
LGE_sweep(Ts, J=J)

J = 1.7
Ts = expspace(-1.2, -1.25, 5)
LGE_sweep(Ts, J=J)

J = 1.8
Ts = expspace(-1.18, -1.23, 5)
LGE_sweep(Ts, J=J)

J = 1.9
Ts = expspace(-1.15, -1.25, 5)
LGE_sweep(Ts, J=J)

J = 2
Ts = expspace(-1.14, -1.19, 5)
LGE_sweep(Ts, J=J)

include("visualize.jl")