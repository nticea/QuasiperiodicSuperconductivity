
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
V0 = 1
V1 = -1.5
LGE_tol = 1e-2
BdG_tol = 1e-12
niter = 500

## Tc using LGE ##

## Superfluid stiffness calculation ##
T = 0.31 # everything is at 0 temperature
Js = collect(0:0.5:4)
cmap = cgrad(:matter, length(Js), categorical=true)
p = plot()
ptime = plot()
pmem = plot()
for (i, J) in enumerate(Js)
    Λs = T * collect(1:1:15)
    λs = []
    times = []
    mems = []
    for Λ in Λs
        res = @timed pairfield_correlation(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, Λ=Λ)
        λ, Δ_LGE = res.value
        push!(λs, λ)
        push!(times, res.time)
        push!(mems, res.bytes)
    end
    res = @timed pairfield_correlation(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)
    λ0, _ = res.value
    t0, m0 = res.time, res.bytes

    plot!(p, Λs, λs, color=cmap[i], label=nothing)
    scatter!(p, Λs, λs, color=cmap[i], label="J=$J")
    hline!(p, [λ0], color=cmap[i], ls=:dashdot, label=nothing)

    plot!(ptime, Λs, times, color=cmap[i], label=nothing)
    scatter!(ptime, Λs, times, color=cmap[i], label=nothing)
    hline!(ptime, [t0], color=cmap[i], ls=:dashdot, label=nothing)

    plot!(pmem, Λs, mems, color=cmap[i], label=nothing)
    scatter!(pmem, Λs, mems, color=cmap[i], label=nothing)
    hline!(pmem, [m0], color=cmap[i], ls=:dashdot, label=nothing)
end
title!(p, "Maximum eigenvalue")
xlabel!(p, "Λ")
ylabel!(p, "λ")

plot!(ptime, ylims=(0, 0.25))
title!(ptime, "Runtime")
xlabel!(ptime, "Λ")
ylabel!(ptime, "Time (s)")

title!(pmem, "Memory usage")
xlabel!(pmem, "Λ")
ylabel!(pmem, "Bytes")

ps = plot(p, ptime, pmem, layout=Plots.grid(1, 3,
        widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1000), plot_title="Keeping states within Λ of εF for V0=$V0, V1=$V1 on $L × $L lattice")