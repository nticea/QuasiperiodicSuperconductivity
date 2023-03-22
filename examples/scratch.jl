## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using Interpolations
using NPZ

L = 9 # the full system is L × L 
t = 1 # hopping 
J = 3
Q = (√5 - 1) / 2
μ = 0
pairing_symmetry = "s-wave"

# temperature & potential 
Js = [0, 1, 2, 3, 4]
V0s = [0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2]
Ts = LinRange(0, 1, 20)

λs = zeros(length(Js), length(V0s), length(Ts))
for (k, J) in enumerate(Js)
    for (i, T) in enumerate(Ts)
        for (j, V0) in enumerate(V0s)
            λ = λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
            @show J, T, V0, λ
            λs[k, j, i] = λ
        end
    end
end

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
title!(p2, "Transition temperature")
xlabel!(p2, "V")
ylabel!(p2, "Tc")


