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
include("../../src/BdG_dwave.jl")
include("../../src/BdG.jl")

## PARAMETERS ##
L = 9 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = nothing#π / 7
ϕx = 0#0.375
ϕy = 0#0.29
V0 = 1
V1 = -1
J = 1
periodic = true
niter = 500
tol = 1e-10
pairing_symmetry = "d-wave"
noise = 0#1e-3 # looks like adding noise is not helping anything 

# T = 0.06
# λ, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")
# @assert 1 == 0

Ts = LinRange(0.008, 0.15, 50)

global λs = []
global Δs = []
global Δijs = []
global Δijs_swave = []
global max_Δs_list = []

for T in Ts
    λ, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")
    Δij, max_Δ = compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, noise=noise)
    if pairing_symmetry == "s-wave"
        Δij_swave, max_Δ_swave = compute_Δ(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, periodic=periodic, niter=niter, tol=tol, noise=noise)
        push!(Δijs_swave, Δij_swave)
    end

    push!(λs, λ)
    push!(Δs, Δ)
    push!(Δijs, Δij)
    push!(max_Δs_list, max_Δ)
end

# convergence data
cmap = cgrad(:matter, 50, categorical=true)
p = plot()
for i in 1:length(Ts)
    plot!(p, max_Δs_list[i], c=cmap[i], label=nothing, yaxis=:log10)
end
xlabel!(p, "Iteration")
ylabel!(p, "Maximum Δ")
title!(p, "Convergence data")

idxlist = sortperm(λs)
knots = λs[idxlist]
Interpolations.deduplicate_knots!(knots, move_knots=true)
interp_linear = linear_interpolation(knots, Ts[idxlist])
Tc = interp_linear(1)
@show Tc

max_Δs = []
for Δij in Δijs
    push!(max_Δs, maximum(Δij))
end

if pairing_symmetry == "s-wave"
    max_Δs_swave = []
    for Δij_swave in Δijs_swave
        push!(max_Δs_swave, maximum(Δij_swave))
    end
end

plot(Ts, max_Δs, label="BdG Δ d-wave", color="blue", yaxis=:log10)
scatter!(Ts, max_Δs, label=nothing, color="blue", yaxis=:log10)
if pairing_symmetry == "s-wave"
    plot!(Ts, max_Δs_swave, label="BdG Δ s-wave", color="red")#, yaxis=:log10)
    scatter!(Ts, max_Δs_swave, label=nothing, color="red")#, yaxis=:log10)
end
vline!([Tc], label="LGE Tc", color="orange")
title!("L=$L, V0=$V0, V1=$V1, J=$J")
xlabel!("T")
ylabel!("Δ")

## PLOT Δ ##
# maxev = Δijs[3]

# function colour_phase(x1::Int, x2::Int, x3::Int; all_evs, numpts::Int=10)
#     cm = palette([:blue, :red], 2 * numpts + 1)
#     val = all_evs[x1, x2, x3]
#     max = maximum(abs.(all_evs))
#     idx = floor(Int, val / max * numpts + numpts + 1)
#     return cm[idx]
# end

# p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
# if λ > 0
#     for x in 1:L
#         for y in 1:L

#             # bonds 
#             plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
#             plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
#             plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
#             plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)

#             # onsite dot 
#             if abs.(maximum(evs[5, x, y])) > 1e-6
#                 scatter!(p, [x], [y], ms=100 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)
#             end

#         end
#     end
# end
# xlabel!(p, "Site (x)")
# ylabel!(p, "Site, (y)")
# title!(p, "T=$T, λ=$(round(λ,digits=2)) \n Δ(J=$J, θ=$θ, V0=$V0, V1=$(round(V1,digits=2)))", fontsize=4)
# plot(p)