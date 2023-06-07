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
L = 17 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
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
        if !already_calculated(df; L=L, J=J, θ=θ, V0=V0, V1=V1, T=T)

            # calculate M at a given (J,T,V0, V1)
            @time λ, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)

            update_results!(df; L=L, λ=λ, J=J, θ=θ, V0=V0, V1=V1, T=T, Δ=Δ)
            CSV.write(savepath, df)
            flush(stdout)
        end
    end
end

J = 0
Ts = expspace(-0.92, -0.8, 5)
LGE_sweep(Ts, J=J)

J = 0.1
Ts = expspace(-0.94, -0.82, 5)
LGE_sweep(Ts, J=J)

J = 0.2
Ts = expspace(-0.96, -0.84, 5)
LGE_sweep(Ts, J=J)

J = 0.3
Ts = expspace(-0.98, -0.86, 5)
LGE_sweep(Ts, J=J)

J = 0.4
Ts = expspace(-1, -0.88, 5)
LGE_sweep(Ts, J=J)

J = 0.5
Ts = expspace(-1.2, -1, 5)
LGE_sweep(Ts, J=J)

J = 0.6
Ts = expspace(-1.4, -1, 5)
LGE_sweep(Ts, J=J)

J = 0.7
Ts = expspace(-1.4, -1, 5)
LGE_sweep(Ts, J=J)

J = 0.8
Ts = expspace(-1.4, -1, 5)
LGE_sweep(Ts, J=J)

J = 0.9
Ts = expspace(-0.94, -0.82, 5)
LGE_sweep(Ts, J=J)

J = 1
Ts = expspace(-0.92, -0.8, 5)
LGE_sweep(Ts, J=J)

J = 1.1
Ts = expspace(-0.9, -0.79, 5)
LGE_sweep(Ts, J=J)

J = 1.2
Ts = expspace(-0.88, -0.78, 5)
LGE_sweep(Ts, J=J)

J = 1.3
Ts = expspace(-0.87, -0.77, 5)
LGE_sweep(Ts, J=J)

J = 1.4
Ts = expspace(-0.86, -0.76, 5)
LGE_sweep(Ts, J=J)

J = 1.5
Ts = expspace(-0.85, -0.75, 5)
LGE_sweep(Ts, J=J)

J = 1.6
Ts = expspace(-0.84, -0.74, 5)
LGE_sweep(Ts, J=J)

J = 1.7
Ts = expspace(-0.83, -0.73, 5)
LGE_sweep(Ts, J=J)

J = 1.8
Ts = expspace(-0.82, -0.72, 5)
LGE_sweep(Ts, J=J)

J = 1.9
Ts = expspace(-0.81, -0.71, 5)
LGE_sweep(Ts, J=J)

J = 2
Ts = expspace(-0.8, -0.7, 5)
LGE_sweep(Ts, J=J)

include("visualize.jl")