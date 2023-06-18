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
θ = nothing #π / 7
ϕx = 0#0.375
ϕy = 0#0.29
V0 = 1
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

L = 17
J = 1
T = 0
λ, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")

@assert 1 == 0

J = 0
Ts = LinRange(0.08, 0.16, 5)
LGE_sweep(Ts, J=J)

J = 0.5
Ts = LinRange(0.07, 0.15, 5)
LGE_sweep(Ts, J=J)

J = 1
Ts = LinRange(0.04, 0.12, 5)
LGE_sweep(Ts, J=J)

# J = 0
# Ts = expspace(-0.62, -0.52, 5)
# LGE_sweep(Ts, J=J)

# J = 0.1
# Ts = expspace(-0.64, -0.54, 5)
# LGE_sweep(Ts, J=J)

# J = 0.2
# Ts = expspace(-0.64, -0.54, 5)
# LGE_sweep(Ts, J=J)

# J = 0.3
# Ts = expspace(-0.62, -0.52, 5)
# LGE_sweep(Ts, J=J)

# J = 0.4
# Ts = expspace(-0.65, -0.55, 5)
# LGE_sweep(Ts, J=J)

# J = 0.5
# Ts = expspace(-0.66, -0.56, 5)
# LGE_sweep(Ts, J=J)

# J = 0.6
# Ts = expspace(-0.67, -0.57, 5)
# LGE_sweep(Ts, J=J)

# J = 0.7
# Ts = expspace(-0.67, -0.57, 5)
# LGE_sweep(Ts, J=J)

# J = 0.8
# Ts = expspace(-0.67, -0.57, 5)
# LGE_sweep(Ts, J=J)

# J = 0.9
# Ts = expspace(-0.66, -0.56, 5)
# LGE_sweep(Ts, J=J)

# J = 1
# Ts = expspace(-0.66, -0.56, 5)
# LGE_sweep(Ts, J=J)

# J = 1.1
# Ts = expspace(-0.65, -0.55, 5)
# LGE_sweep(Ts, J=J)

# J = 1.2
# Ts = expspace(-0.64, -0.54, 5)
# LGE_sweep(Ts, J=J)

# J = 1.3
# Ts = expspace(-0.64, -0.54, 5)
# LGE_sweep(Ts, J=J)

# J = 1.4
# Ts = expspace(-0.64, -0.54, 5)
# LGE_sweep(Ts, J=J)

# J = 1.5
# Ts = expspace(-0.63, -0.53, 5)
# LGE_sweep(Ts, J=J)

# J = 1.6
# Ts = expspace(-0.63, -0.53, 5)
# LGE_sweep(Ts, J=J)

# J = 1.7
# Ts = expspace(-0.63, -0.53, 5)
# LGE_sweep(Ts, J=J)

# J = 1.8
# Ts = expspace(-0.62, -0.52, 5)
# LGE_sweep(Ts, J=J)

# J = 1.9
# Ts = expspace(-0.61, -0.51, 5)
# LGE_sweep(Ts, J=J)

# J = 2
# Ts = expspace(-0.61, -0.51, 5)
# LGE_sweep(Ts, J=J)

# include("visualize.jl")