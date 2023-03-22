## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using Interpolations
using NPZ

L = 5 # the full system is L × L 
t = 1 # hopping 
J = 3
Q = (1 + √5) / 2
μ = 0
pairing_symmetry = "s-wave"

# temperature & potential 
V0s = [0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]
Ts = vcat(LinRange(0, 0.49, 50), LinRange(0.49, 1, 50))

λs1 = npzread("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/results.npy")
λs2 = npzread("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/results2.npy")
λs3 = npzread("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/results3.npy")

λs = NaN * zeros(length(V0s), length(Ts))
λs[1:14, 1:50] .= λs1
λs[15:end, 1:50] .= λs2
λs[15:end, 51:end] .= λs3

# visualize everything 
global p = plot()
cmap = cgrad(:acton, length(V0s), categorical=true)
for (j, V0) in enumerate(V0s)
    global p = plot!(p, Ts, λs[j, :], label=nothing, c=cmap[j])
    global p = scatter!(p, Ts, λs[j, :], label="V=$(V0)", c=cmap[j])
end
plot(p)

# do interpolations 
Tcs = []
for j in 1:14
    knots = reverse(λs[j, 1:50])
    Interpolations.deduplicate_knots!(knots)
    try
        interp_linear = linear_interpolation(knots, reverse(Ts[1:50]))
        push!(Tcs, interp_linear(1))
    catch e
        push!(Tcs, NaN)
    end
end
for j in 15:length(V0s)
    knots = reverse(λs[j, :])
    Interpolations.deduplicate_knots!(knots)
    try
        interp_linear = linear_interpolation(knots, reverse(Ts))
        push!(Tcs, interp_linear(1))
    catch e
        push!(Tcs, NaN)
    end
end
p2 = plot(V0s, Tcs, xaxis=:log10, yaxis=:log10, c="blue", label=nothing)
scatter!(p2, V0s, Tcs, xaxis=:log10, yaxis=:log10, c="blue", label="J=3")
title!(p2, "Transition temperature")
xlabel!(p2, "V")
ylabel!(p2, "Tc")

