
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")

## PARAMETERS ##
L = 17 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
periodic = true
ϕx = 0
ϕy = 0
LGE_tol = 1e-2
BdG_tol = 1e-12
niter = 500

J = 0

symmetry = "d-wave"
Ts = expspace(-5, 1, 20)
Ls = [23]

if symmetry == "s-wave"
    χs = zeros(length(Ts), length(Ls))
else
    χs = zeros(length(Ts), length(Ls), 3, 3)
end

for (i, T) in enumerate(Ts)
    for (j, L) in enumerate(Ls)
        χ = @time pairfield_susceptibility(T, symmetry, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic)
        if symmetry == "s-wave"
            χs[i, j] = χ
        else
            χs[i, j, :, :] = χ
        end
    end
end

p = plot()
cmap = ["blue"]#cgrad(:matter, length(Ls), categorical=true)
for (j, L) in enumerate(Ls)
    if symmetry == "s-wave"
        plot!(Ts, χs[:, j], color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(Ts, χs[:, j], color=cmap[j], label="L=$L", xaxis=:log10)
    else
        plot!(Ts, χs[:, j, 1, 1], color=cmap[j], ls=:dash, label=nothing, xaxis=:log10)
        scatter!(Ts, χs[:, j, 1, 1], color=cmap[j], label="x̂, L=$L", xaxis=:log10)
        plot!(Ts, χs[:, j, 2, 2], color=cmap[j], ls=:dashdot, label=nothing, xaxis=:log10)
        scatter!(Ts, χs[:, j, 2, 2], color=cmap[j], label="ŷ, L=$L", xaxis=:log10)
        plot!(Ts, χs[:, j, 3, 3], color="red", label=nothing, xaxis=:log10)
        scatter!(Ts, χs[:, j, 3, 3], color="red", label="on-site, L=$L", xaxis=:log10)
    end
end
title!("Susceptibility for J=$J")
xlabel!("T")
ylabel!("χ")

# plot(Ts, χs, color="blue", label=nothing, xaxis=:log10)
# scatter!(Ts, χs, color="blue", label="$symmetry", xaxis=:log10)
# title!("Susceptibility for $L × $L system for J=$J")
# xlabel!("T")
# ylabel!("χ")