
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")

## PARAMETERS ##
L = 17
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
periodic = true
ϕx = 0
ϕy = 0
LGE_tol = 1e-2
BdG_tol = 1e-12
niter = 500

symmetry = "s-wave"
J = 0
T = 1

χ_true = @time pairfield_susceptibility(T, symmetry, L=L, t=t, J=J,
    Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic, Λ=nothing)

# χ_fast = @time uniform_susceptibility(T, symmetry, L=L, t=t, J=J,
#     Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic, Λ=nothing)

Ts = expspace(-3, 1, 20)
Ls = [17, 23, 31, 35]
χs = zeros(length(Ls), length(Ts))
for (i, T) in enumerate(Ts)
    for (j, L) in enumerate(Ls)
        χ = @time uniform_susceptibility(T, symmetry, L=L, t=t, J=J,
            Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic, Λ=nothing)
        if symmetry == "d-wave"
            # make the d-wave components 
            xx = χ[1, 1]
            yy = χ[2, 2]
            xy = χ[1, 2]
            yx = χ[2, 1]
            χ = xx + yy - xy - yx
        end
        χs[j, i] = χ / (L * L)
    end
end

p = plot()
c = palette(:tab10)
for (j, L) in enumerate(Ls)
    plot!(p, Ts, χs[j, :], c=c[j], label=nothing, xaxis=:log10)
    scatter!(p, Ts, χs[j, :], c=c[j], label="L=$L", xaxis=:log10)
end
title!("$symmetry susceptibility")
xlabel!("T")
ylabel!("χ")