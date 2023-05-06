## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using LinearAlgebra, Arpack, Plots
using Profile
include("../src/model.jl")
include("../src/meanfield.jl")
# include("../src/BdG_dwave.jl")
# include("../src/BdG.jl")

## PARAMETERS ##
L = 20 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = nothing
J = 3
T = 0.24
V0 = 1
V1 = -1.5
periodic = true

@time M = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="s-wave")
@time M = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")

@time decomp, _ = partialschur(Hermitian(M), nev=1, tol=1e-6, which=LM())
maxev = decomp.Q
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
        plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)

        # onsite dot 
        # scatter!(p, [y], [-x], ms=100 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)
        # plot!(p, [y, y], [-x, -x + 1], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
        # plot!(p, [y, y + 1], [-x, -x], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
        # plot!(p, [y, y], [-x, -x - 1], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
        # plot!(p, [y, y - 1], [-x, -x], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)
    end
end
plot(p)
xlabel!(p, "Site (x)")
ylabel!(p, "Site, (y)")
title!("Δ(J=$J, V0=$V0, V1=$V1)")

# , c=colour_phase(evs[5, x, y])


# visualization
# scatter plot each of the points (the last block). Size is determined by magnitude, and colour by phase
# draw lines between the points like: plot([x1 y1],[x2 y2], linewidth=2,legend=:false)
# control the line width using the value of the eigenvector 


# Δs = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")
# heatmap(Δs[5, :, :])

# Ts = Ts = expspace(-4, 0, 20)
# V1s = [-1, -0.5, 0, 0.5, 1]
# # Δs = zeros(length(Ts), length(V1s))
# Tcs = zeros(length(V1s))
# for (j, V1) in enumerate(V1s)
#     λs = []
#     for (i, T) in enumerate(Ts)
#         λ = λmax(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")
#         push!(λs, λ)
#         @show λ, V1, λ
#         # Δ = compute_Δ_dwave(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, tol=tol, periodic=periodic)
#         # @show T, V1, Δ
#         # Δs[i, j] = Δ
#     end
#     knots_λs = reverse(λs)
#     Interpolations.deduplicate_knots!(knots_λs, move_knots=true)
#     interp_linear = linear_interpolation(knots_λs, reverse(Ts))
#     Tc = interp_linear(1)
#     push!(Tcs, Tc)
# end

# p = plot()
# cs = palette([:purple, :black], length(V1s))
# for (i, V1) in enumerate(V1s)
#     global p = plot(p, Ts, Δs[:, i], c=cs[i], label=nothing, xscale=:log10)
#     global p = scatter(p, Ts, Δs[:, i], c=cs[i], label="V1=$V1", xscale=:log10)
#     #vline!(p, [Tcs[i]], c=cs[i], label="V1=$V1")
# end
# xlabel!(p, "T1")
# ylabel!(p, "Δ")


# λs, Δs = [], []
# for T in Ts
#     λ = λmax(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")
#     @show λ
#     Δ = compute_Δ_dwave(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, tol=tol, periodic=periodic)
#     @show Δ

#     push!(λs, λ)
#     push!(Δs, Δ)
# end

# # plot(Ts, λs, label="λ")

# knots_λs = reverse(λs)
# Interpolations.deduplicate_knots!(knots_λs, move_knots=true)
# interp_linear = linear_interpolation(knots_λs, reverse(Ts))
# Tc = interp_linear(1)

# plot(Ts, Δs, c="blue", label=nothing, xscale=:log10)
# scatter!(Ts, Δs, label="Δ", c="blue", xscale=:log10)
# vline!([Tc], label="Tc from LGE", xscale=:log10)
# title!("Δ from BdG equations, L=$L J=$J")
# xlabel!("T")
# ylabel!("Δ")

# Δij, conv, max_Δ = compute_Δ_dwave(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, periodic=periodic)
# heatmap(reverse(M, dims=2), c=:bwr, clims=(-maximum(abs.(M)), maximum(abs.(M))))
# @show maximum(abs.(Δij))

# Δi, conv, max_Δ = compute_Δ(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, tol=tol, periodic=periodic)
# plot(max_Δ, label=nothing, c="blue")
# scatter!(max_Δ, label="Δ", c="blue", ms=0.8)
# xlabel!("Iteration")
# ylabel!("Δ")
# title!("Convergence for L=$L, J=$J")