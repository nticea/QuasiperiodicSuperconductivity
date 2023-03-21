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
V0s = [0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2]
Ts = LinRange(0, 0.49, 50)

λs = zeros(length(V0s), length(Ts))
for (i, T) in enumerate(Ts)
    for (j, V0) in enumerate(V0s)
        λ = λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
        @show T, λ
        λs[j, i] = λ
    end
end

# save the results 
npzwrite("results.npy", λs)

# visualize everything 
global p = plot()
cmap = palette(:tab10)
for (j, V0) in enumerate(V0s)
    global p = plot!(p, Ts, λs[j, :], label=nothing, c=cmap[j])
    global p = scatter!(p, Ts, λs[j, :], label="V=$(V0)", c=cmap[j])
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
plot(V0s, Tcs, xaxis=:log10, yaxis=:log10)


