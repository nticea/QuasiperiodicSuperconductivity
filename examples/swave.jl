## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using Interpolations
using NPZ

L = 5 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
pairing_symmetry = "s-wave"

# other parameters
Js = [0, 1, 2, 3, 4]
V0s = collect(0.3:0.015:2)
Ts = LinRange(0, 3.5, 40)

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

#save the results 
savepath = joinpath(@__DIR__, "manyJ_results.npy")
npzwrite(savepath, λs)

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


