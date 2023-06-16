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
ϕx = 0.375
ϕy = 0.29
V0 = -2.3
V1 = 0
periodic = true
pairing_symmetry = "s-wave"

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_swave_results.csv")

# load in the dataframe, if it exists. If not, make a new one
df = load_dataframe(savepath)

function LGE_sweep(Ts; J::Real)
    ## RUNNING THE CODE ## 
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))

    for i in iter # iterate through all temperatures
        T = Ts[i]
        # check whether this particular (J,T,V0) combo has been already computed 
        if !already_calculated(df; L=L, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, T=T)

            # calculate M at a given (J,T,V0, V1)
            @time λ, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)

            update_results!(df; L=L, λ=λ, J=J, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, T=T, Δ=Δ)
            CSV.write(savepath, df)
            flush(stdout)
        end
    end
end

J = 0
Ts = expspace(-0.62, -0.52, 5)
LGE_sweep(Ts, J=J)

J = 0.1
Ts = expspace(-0.6, -0.5, 5)
LGE_sweep(Ts, J=J)

J = 0.2
Ts = expspace(-0.6, -0.5, 5)
LGE_sweep(Ts, J=J)

J = 0.3
Ts = expspace(-0.6, -0.5, 5)
LGE_sweep(Ts, J=J)

J = 0.4
Ts = expspace(-0.55, -0.45, 5)
LGE_sweep(Ts, J=J)

J = 0.5
Ts = expspace(-0.5, -0.4, 5)
LGE_sweep(Ts, J=J)

J = 0.6
Ts = expspace(-0.5, -0.4, 5)
LGE_sweep(Ts, J=J)

J = 0.7
Ts = expspace(-0.5, -0.4, 5)
LGE_sweep(Ts, J=J)

J = 0.8
Ts = expspace(-0.5, -0.4, 5)
LGE_sweep(Ts, J=J)

J = 0.9
Ts = expspace(-0.5, -0.4, 5)
LGE_sweep(Ts, J=J)

J = 1
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.1
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.2
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.3
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.4
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.5
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.6
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.7
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.8
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 1.9
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

J = 2
Ts = expspace(-0.45, -0.35, 5)
LGE_sweep(Ts, J=J)

include("visualize.jl")
