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
ϕx = 0
ϕy = 0
V0 = 1
V1 = -1.5
periodic = true

Js = collect(0:0.1:2)


# α = 1.5
# V0s = LinRange(0.1, 2, 10)
# Js = collect(0:0.1:2)

# Tcs_dwave = zeros(length(Js), length(V0s))
# for (i, J) in enumerate(Js)
#     for (k, V0) in enumerate(V0s)
#         V1 = -V0 / α
#         @show (J, V0, V1)
#         Tcs_dwave[i, k] = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry, tol=1e-3)
#     end
# end

# p = plot(size=(750, 500))
# cmap = cgrad(:matter, length(V0s), categorical=true)
# for (k, V0) in enumerate(V0s)
#     V1 = -V0 / α
#     plot!(p, Js, Tcs_dwave[:, k], label=nothing, color=cmap[k])
#     scatter!(p, Js, Tcs_dwave[:, k], label="V0=$(round(V0, digits=2)), V1=$(round(V1, digits=2))", color=cmap[k])
# end
# title!(p, "Transition temperature for |V0/V1|=$α on $L×$L lattice", fontsize=12)
# xlabel!(p, "J")
# ylabel!(p, "Tc")
# plot!(p, legend=:top)

# α = 1.5
# V0s = [1.5]
# Js = collect(0:0.1:2)
# Q = (√5 - 1) / 2
# V0_swave = -1.75

# Tcs_dwave = zeros(length(Js), length(V0s))
# Tcs_swave = zeros(length(Js), length(V0s))
# for (i, J) in enumerate(Js)
#     for (k, V0) in enumerate(V0s)
#         V1 = -V0 / α
#         @show (J, V0, V1)
#         Tcs_dwave[i, k] = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave", tol=1e-3)
#         Tcs_swave[i, k] = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0_swave, V1=V1_swave, periodic=periodic, symmetry="s-wave", tol=1e-3)
#     end
# end

# p = plot(size=(750, 500))
# for (k, V0) in enumerate(V0s)
#     plot!(p, Js, Tcs_dwave[:, k], label=nothing, color="blue")
#     scatter!(p, Js, Tcs_dwave[:, k], label="dwave", color="blue")
#     plot!(p, Js, Tcs_swave[:, k], label=nothing, color="red")
#     scatter!(p, Js, Tcs_swave[:, k], label="swave", color="red")
# end
# title!(p, "Transition temperature for |V0/V1|=$α on $L×$L lattice", fontsize=12)
# xlabel!(p, "J")
# ylabel!(p, "Tc")
# plot!(p, legend=:top)