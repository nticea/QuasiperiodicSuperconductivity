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
Q = (√5 - 1) / 2
μ = 0
pairing_symmetry = "s-wave"

# temperature & potential 
V0s = collect(0.3:0.015:7)
Ts = LinRange(0, 3.5, 40)

λs = npzread("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/results_full.npy")

# visualize everything 
global p = plot()
cmap = cgrad(:acton, length(V0s), categorical=true)
for (j, V0) in enumerate(V0s)
    global p = plot!(p, Ts, λs[j, :], label=nothing, c=cmap[j])
    #global p = scatter!(p, Ts, λs[j, :], label="V=$(V0)", c=cmap[j])
end
plot(p)

# do interpolations 
Tcs = []
for j in 1:length(V0s)
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

