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
J = 0
V1 = 0.5
V0 = 1
periodic = true
pairing_symmetry = "d-wave"

T = 0.04
M = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)

# perform the decomposition 
decomp, _ = partialschur(Hermitian(M), nev=1, tol=1e-6, which=LM())
# extract the maximum eigenvector/value pair 
maxev = decomp.Q[:, 1]
λ = decomp.R[1]
@show λ

evs = zeros(5, L, L)
for (n, i) in enumerate(1:(L*L):(5*L*L))
    evi = maxev[i:(i+L*L-1)]
    evs[n, :, :] = reshape(evi, L, L)
end

function colour_phase(x1::Int, x2::Int, x3::Int; all_evs, numpts::Int=10)
    cm = palette([:blue, :red], 2 * numpts + 1)
    val = all_evs[x1, x2, x3]
    max = maximum(abs.(all_evs))
    idx = floor(Int, val / max * numpts + numpts + 1)
    return cm[idx]
end

p = plot(xlims=(0, L + 1), ylims=(0, L + 1))
# p = plot(xlims=(0, L + 1), ylims=(-L - 1, 0))
for x in 1:L
    for y in 1:L
        # onsite dot 
        scatter!(p, [x], [y], ms=100 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)

        # bonds 
        plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)
    end
end
xlabel!(p, "Site (x)")
ylabel!(p, "Site, (y)")
title!(p, "T=$T, λ=$(round(λ,digits=2)): Δ(J=$J, θ=$θ, V0=$V0, V1=$(round(V1,digits=2)))", fontsize=6)
return p