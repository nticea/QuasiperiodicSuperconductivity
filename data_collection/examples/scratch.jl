
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
L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
ϕx = 0
ϕy = 0
V0 = 1
V1 = -1.5
periodic = true
niter = 1000
tol = 1e-15
T = 0

χs_list_LGE = []
χs_list_BdG = []
χs_list_BdG_primed = []

function primed_Δ(Δ)

    if ndims(Δ) == 1
        Δ = spatial_profile(Δ, L=L)
    end

    Δrot = zeros(size(Δ)...)
    for x in 1:L
        for y in 1:L
            for n in 1:5
                x̃ = -y + L + 1
                ỹ = x
                Δrot[n, x̃, ỹ] = Δ[n, x, y]
            end
        end
    end

    # now swap bonds
    Δrot0 = copy(Δrot)
    Δrot[1, :, :] = Δrot0[2, :, :]
    Δrot[2, :, :] = Δrot0[3, :, :]
    Δrot[3, :, :] = Δrot0[4, :, :]
    Δrot[4, :, :] = Δrot0[1, :, :]

    # take the difference 
    Δ_diff = Δ - Δrot

    # normalize 
    Δ_diff .*= norm(Δ) / norm(Δ_diff)

    return Δ_diff
end

Js = collect(0:0.1:2)
for J in Js
    λ, Δ_LGE = pairfield_correlation(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")
    push!(χs_list_LGE, symmetry_character(Δ_LGE, L=L))
    Δ_BdG, hist = compute_Δ_dwave(T; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, Δ_init=Δ_LGE)
    push!(χs_list_BdG, symmetry_character(Δ_BdG, L=L))
    primed_Δ(Δ_BdG)
end

plot(Js, χs_list_LGE, color="red", label=nothing)
scatter!(Js, χs_list_LGE, color="red", label="LGE soln")
plot!(Js, χs_list_BdG, color="blue", label=nothing)
scatter!(Js, χs_list_BdG, color="blue", label="BdG soln")
title!("Symmetry parameter for \n Δ(V0=$V0, V1=$V1, θ=$(θ_to_π(θ)),ϕx=$ϕx,ϕy=$ϕy) \n $L × $L lattice at T=$(round(T,digits=2))")
xlabel!("J")