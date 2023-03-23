## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using Interpolations
using NPZ

L = 9 # the full system is L × L 

# temperature & potential 
Js = [0, 1, 2, 3, 4]
V0s = [0.37, 0.38, 0.39, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2]
Ts = LinRange(0, 1, 20)

λs = npzread("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/9Nx9Ny_results.npy")

# do interpolations 
Tcs = zeros(length(Js), length(V0s))
for k in 1:length(Js)
    for j in 1:length(V0s)
        knots = reverse(λs[k, j, :])
        Interpolations.deduplicate_knots!(knots)
        try
            interp_linear = linear_interpolation(knots, reverse(Ts))
            Tcs[k, j] = interp_linear(1)
        catch e
            Tcs[k, j] = NaN
        end
    end
end

p2 = plot()
cmap = cgrad(:Set1_9, length(V0s), categorical=true)
for (k, J) in enumerate(Js)
    plot!(p2, V0s, Tcs[k, :], xaxis=:log10, yaxis=:log10, c=cmap[k], label=nothing)
    scatter!(p2, V0s, Tcs[k, :], xaxis=:log10, yaxis=:log10, c=cmap[k], label="J=$(J)")
end
title!(p2, "Transition temperature for $(L)x$(L) square lattice")
xlabel!(p2, "V")
ylabel!(p2, "Tc")

