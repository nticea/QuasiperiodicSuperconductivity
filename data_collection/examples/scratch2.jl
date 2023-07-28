
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
L = 17
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
periodic = true
ϕx = 0
ϕy = 0
J = 0
Λ = 1

Ts = expspace(-5, 1, 10)
Ls = [55]
Λs = [0.5, 1, 2, 3, nothing]
χs_dwave = zeros(length(Λs), length(Ts))
χs_swave = zeros(length(Λs), length(Ts))

# χ = @time uniform_susceptibility(1e-2, "d-wave", L=67, t=t, J=J,
#     Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic, Λ=1)

# @assert 1 == 0

for (i, T) in enumerate(Ts)
    for (j, Λ) in enumerate(Λs)
        χ = @time uniform_susceptibility(T, L=L, t=t, J=J,
            Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic, Λ=Λ)
        # make the d-wave components 
        xx = χ[1, 1]
        yy = χ[2, 2]
        xy = χ[1, 2]
        yx = χ[2, 1]

        χs_dwave[j, i] = xx + yy - xy - yx
        χs_swave[j, i] = χ[3, 3]
    end
end

pswave = plot()
c = palette(:tab10)
for (j, Λ) in enumerate(Λs)
    plot!(pswave, Ts, χs_swave[j, :], c=c[j], label=nothing, xaxis=:log10)
    scatter!(pswave, Ts, χs_swave[j, :], c=c[j], label="Λ=$Λ", xaxis=:log10)
end
title!("s-wave susceptibility")
xlabel!("T")
ylabel!("χ")

pdwave = plot()
c = palette(:tab10)
for (j, Λ) in enumerate(Λs)
    plot!(pdwave, Ts, χs_dwave[j, :], c=c[j], label=nothing, xaxis=:log10)
    scatter!(pdwave, Ts, χs_dwave[j, :], c=c[j], label="Λ=$Λ", xaxis=:log10)
end
title!("d-wave susceptibility")
xlabel!("T")
ylabel!("χ")

p = plot(pswave, pdwave, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 600), plot_title="LGE χ for $Λ cutoff at J=$J")