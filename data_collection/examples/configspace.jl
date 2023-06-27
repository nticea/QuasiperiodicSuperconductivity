## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/model.jl")
include("../../src/meanfield.jl")
include("../../src/BdG_dwave.jl")

## PARAMETERS ##
L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0
ϕy = 0
V0 = 1.5
V1 = -1
J = 1.5
periodic = true
pairing_symmetry = "d-wave"
niter = 500
tol = 1e-12

plot_potential(L=L, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy)

T = 0
λ, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)
Δ_BdG, hist = compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, Δ_init=Δ_LGE)

Δ = spatial_profile(Δ_LGE, L=L)
p0 = plot_spatial_profile(Δ, L=L, title="Real space")
p1 = plot_in_config_space(Δ[5, :, :], L=L, Q=Q, θ=θ, title="On-site")
p2 = plot_in_config_space(Δ[1, :, :], L=L, Q=Q, θ=θ, title="-x̂ bond")
p3 = plot_in_config_space(Δ[2, :, :], L=L, Q=Q, θ=θ, title="+ŷ bond")
p4 = plot_in_config_space(Δ[3, :, :], L=L, Q=Q, θ=θ, title="x̂ bond")
p5 = plot_in_config_space(Δ[4, :, :], L=L, Q=Q, θ=θ, title="-ŷ bond")
p = plot(p1, p2, p3, p4, p5, p0, layout=Plots.grid(2, 3,
        widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1000), aspect_ratio=:equal, plot_title=" LGE Δ(J=$J, V0=$V0, V1=$V1, θ=$(θ_to_π(θ))) for $L × $L lattice at T=$(round(T,digits=2))")

Δ = spatial_profile(Δ_BdG, L=L)
p0 = plot_spatial_profile(Δ, L=L, title="Real space")
p1 = plot_in_config_space(Δ[5, :, :], L=L, Q=Q, θ=θ, title="On-site")
p2 = plot_in_config_space(Δ[1, :, :], L=L, Q=Q, θ=θ, title="-x̂ bond")
p3 = plot_in_config_space(Δ[2, :, :], L=L, Q=Q, θ=θ, title="+ŷ bond")
p4 = plot_in_config_space(Δ[3, :, :], L=L, Q=Q, θ=θ, title="x̂ bond")
p5 = plot_in_config_space(Δ[4, :, :], L=L, Q=Q, θ=θ, title="-ŷ bond")
p = plot(p1, p2, p3, p4, p5, p0, layout=Plots.grid(2, 3,
        widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1000), aspect_ratio=:equal, plot_title=" BdG Δ(J=$J, V0=$V0, V1=$V1, θ=$(θ_to_π(θ)),ϕx=$ϕx,ϕy=$ϕy) for $L × $L lattice at T=$(round(T,digits=2))")

