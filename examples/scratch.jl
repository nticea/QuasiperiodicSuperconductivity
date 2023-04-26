## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
include("../src/BdG_dwave.jl")

using Plots
using ProgressBars
using CSV
using DataFrames

## PARAMETERS ##
L = 20 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = nothing
J = 2
T = 0
V0 = 1.5
V1 = 0
periodic = true
tol = 1e-4

# run the BdG code 
Δ, conv, max_Δ = compute_Δ_dwave(T, L=L, t=t, J=J, Q=Q, μ=μ, periodic=periodic, V0=V0, V1=V1, θ=θ, tol=tol)

heatmap(reverse(Δ, dims=2), cmap=:bwr, clims=(-maximum(abs.(Δ)), maximum(abs.(Δ))))
title!("Δ_ij for $(L)x$(L) system")

fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, periodic=periodic))
plot(max_Δ, label="Δ_ij (max)")
xlabel!("Iteration")
ylabel!("Maximum Δ_ij")
title!("Convergence of maximum Δ_ij")
hline!([fsgap], label="Finite size gap")