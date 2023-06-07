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
L = 9 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
J = 0.5
V1 = -1
V0 = 1.5
periodic = true
pairing_symmetry = "d-wave"

T = 0.12
λ_h, Δ_h = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)
λ_nh, Δ_nh = pairfield_correlation_nonhermitian(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)

maxev = Δ_nh
λ = λ_nh
θ = θ_to_π(θ)

if pairing_symmetry == "d-wave"
    evs = zeros(5, L, L)
    for (n, i) in enumerate(1:(L*L):(5*L*L))
        evi = maxev[i:(i+L*L-1)]
        evs[n, :, :] = reshape(evi, L, L)
    end
elseif pairing_symmetry == "s-wave"
    evs = reshape(maxev, L, L)
end

function colour_phase(x1::Int, x2::Int, x3::Int; all_evs, numpts::Int=10)
    cm = palette([:blue, :red], 2 * numpts + 1)
    val = all_evs[x1, x2, x3]
    max = maximum(abs.(all_evs))
    idx = floor(Int, val / max * numpts + numpts + 1)
    return cm[idx]
end

p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
if λ > 0
    for x in 1:L
        for y in 1:L

            # bonds 
            plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)

            # onsite dot 
            if abs.(maximum(evs[5, x, y])) > 1e-6
                scatter!(p, [x], [y], ms=100 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)
            end

        end
    end
end
xlabel!(p, "Site (x)")
ylabel!(p, "Site, (y)")
title!(p, "T=$T, λ=$(round(λ,digits=2)) \n Δ(J=$J, θ=$θ, V0=$V0, V1=$(round(V1,digits=2)))", fontsize=4)
plot(p)