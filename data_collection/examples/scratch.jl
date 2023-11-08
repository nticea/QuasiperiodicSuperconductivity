## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")
include("../../src/stiffness.jl")

# Parameters 
t = 1
L = 7
Q = (√5 - 1) / 2
μ = 0.75
θ = π / 7
V0 = 3
V1 = -2
ϕx = 0
ϕy = 0
ϕz = 0
periodic = true
disorder = false
ndims = 3
T = 0
J = 2
slice = 1

# simulation parameters 
niter = 500
BdG_tol = 1e-15
LGE_tol = 1e-2

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
# λ, Δ_LGE = pairfield_correlation(m, T=T)
Tc, λ, Δ_LGE = LGE_find_Tc(m, npts=5, tol=LGE_tol)
K, Π, Δ_BdG = superfluid_stiffness_finiteT(m, T=0, tol=BdG_tol, niter=niter, Δ_init=Δ_LGE)

# λ, Δ_LGE = @time pairfield_correlation(m, T=T)

# @assert 1 == 0

# E, U = diagonalize_hamiltonian(m)

# # s-wave case is faster 
# if V1 == 0
#     M = swave(m, T, E=E, U=U)
#     λ, Δ = calculate_λ_Δ(M)
# else
#     M = dwave(m, T, E=E, U=U)
#     λ, Δ = calculate_λ_Δ(M)
# end
# @assert abs(imag(λ)) < 1e-6
# λ = real(λ)

# @assert 1 == 0

Δ = Δ_LGE
Δ = spatial_profile(m, Δ=real.(Δ))
evs = real.(Δ)[:, :, :, slice]

p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false, aspect_ratio=true)
for x in 1:L
    for y in 1:L

        # bonds 
        plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[5, x, y]), alpha=10 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)

        # onsite dot 
        scatter!(p, [x], [y], ms=100 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
    end
end

xlabel!(p, "Site (x)")
ylabel!(p, "Site, (y)")
title!(p, "J=$(m.J), V0=$(m.V0), V1=$(m.V1), θ=$(θ_to_π(m.θ)) on $(m.L)×$(m.L)×$(m.L) lattice at z=$slice")