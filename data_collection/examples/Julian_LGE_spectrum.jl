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
Q0 = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0.375
ϕy = 0.29
periodic = true

T = 0.13332625935631917
V0_dwave = 1.5
V1_dwave = -1
V0_swave = -1.745
V1_swave = 0

Js = collect(0:0.1:2)
max_λ_dwave_Q0 = []
max_λ_swave_Q0 = []
for J in Js
    λ_dwave, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q0, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0_dwave, V1=V1_dwave, periodic=periodic, symmetry="d-wave")
    λ_swave, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q0, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0_swave, V1=V1_swave, periodic=periodic, symmetry="s-wave")

    push!(max_λ_dwave_Q0, λ_dwave)
    push!(max_λ_swave_Q0, λ_swave)
end

max_λ_dwave_4Q0 = []
max_λ_swave_4Q0 = []
for J in Js
    λ_dwave, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=4 * Q0, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0_dwave, V1=V1_dwave, periodic=periodic, symmetry="d-wave")
    λ_swave, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=4 * Q0, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0_swave, V1=V1_swave, periodic=periodic, symmetry="s-wave")

    push!(max_λ_dwave_4Q0, λ_dwave)
    push!(max_λ_swave_4Q0, λ_swave)
end

# plot(Js, max_λ_dwave_Q0, c="blue", ls=:dash, label="d-wave, Q0")
# plot!(Js, max_λ_swave_Q0, c="red", ls=:dash, label="s-wave, Q0")
# plot!(Js, max_λ_dwave_4Q0, c="blue", ls=:dot, label="d-wave, 4 × Q0")
# plot!(Js, max_λ_swave_4Q0, c="red", ls=:dot, label="s-wave, 4 × Q0")
# xlabel!("J")
# title!("Maximum eigenvalue")

plot(Js, max_λ_dwave_Q0 ./ max_λ_swave_Q0, c="black", ls=:dash, label="Q0")
plot!(Js, max_λ_dwave_4Q0 ./ max_λ_swave_4Q0, c="black", ls=:dot, label="4 × Q0")
xlabel!("J")
title!("d-wave/s-wave ratio of max eigenvalue")