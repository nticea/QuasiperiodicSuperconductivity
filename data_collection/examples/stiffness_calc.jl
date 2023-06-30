## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/stiffness.jl")

## PARAMETERS ##
L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0#0.123
ϕy = 0#0.456
V0 = -1.5
V1 = 0
J = 0
periodic = true
niter = 500
tol = 1e-12
noise = 0
npts = 5

@error "Please check factors of 2 in the kinetic, Aq, and Dq terms!"

Ts = expspace(-1, -9, npts)
Ds = zeros(npts, 4)
# perform extrapolation q → 0
for (i, T) in enumerate(Ts)
    print(i, "-")
    K, Π = superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise)
    @show K, Π
    @show -K + Π
    Ds[i, :] = -K + Π
end

# perform extrapolation T → 0
Ds_extrapolated = []
for x in 1:4
    model = Polynomials.fit(Ts, Ds[:, x], npts)
    push!(Ds_extrapolated, model(0))
end

# T = 1e-9
# K, Π = superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise)