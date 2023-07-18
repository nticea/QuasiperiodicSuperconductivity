
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/meanfield.jl")

## PARAMETERS ##
L = 17 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
periodic = true
ϕx = 0.123
ϕy = 0.987
V0 = 1
V1 = -1.5
J = 2.5

Tc, λ, Δ = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, tol=1e-2, npts=5)
p0 = plot_spatial_profile(Δ, L=L, title="Real space LGE Δ")
p1 = plot_potential(L=L, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy)

p = plot(p1, p0, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1200, 700), aspect_ratio=:equal, plot_title="V0=$V0, V1=$V1 on $L × $L lattice at Tc=$(round(Tc,digits=2))")

#savefig(p, "potential.pdf")