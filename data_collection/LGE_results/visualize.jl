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

loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/LGE_results/31Nx31Ny_dwave_results.csv"
idx = 6

df = load_dataframe(loadpath)
L = df.L[idx]
J = df.J[idx]
V0 = df.V0[idx]
V1 = df.V1[idx]
θ = θ_to_π(df.θ[idx])

if length(df.Δ[1]) == 5 * L^2
    symmetry = "d-wave"
elseif length(df.Δ[1]) == L^2
    symmetry = "s-wave"
else
    @error "Symmetry not recognized"
    @assert 1 == 0
end

maxev = df.Δ[idx]
λ = df.λ[idx]
if symmetry == "d-wave"
    evs = zeros(5, L, L)
    for (n, i) in enumerate(1:(L*L):(5*L*L))
        evi = maxev[i:(i+L*L-1)]
        evs[n, :, :] = reshape(evi, L, L)
    end
elseif symmetry == "s-wave"
    evs = reshape(maxev, L, L)
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
title!(p, "$symmetry, λ=$(round(λ,digits=3)): Δ(J=$J, θ=$θ, V0=$V0, V1=$V1)")
plot(p)
