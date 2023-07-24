
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
J = 1.1
Tc, λ, Δ_LGE = @time LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, tol=LGE_tol, npts=5)

function to_5N_LGE_Δ(maxev; L)
    @assert length(maxev) == 3 * L * L
    evs = zeros(3, L, L)
    for (n, i) in enumerate(1:(L*L):(3*L*L))
        evi = maxev[i:(i+L*L-1)]
        evs[n, :, :] = reshape(evi, L, L)
    end

    # now put in the missing blocks
    Δ1 = evs[1, :, :]
    Δ2 = evs[2, :, :]
    Δ3 = circshift(Δ1, (-1, 0)) # mapping between +x̂ and -x̂
    Δ4 = circshift(Δ2, (0, 1)) # mapping between +ŷ and -ŷ
    new_evs = zeros(5, L, L)
    new_evs[1, :, :] = Δ1
    new_evs[2, :, :] = Δ2
    new_evs[3, :, :] = Δ3
    new_evs[4, :, :] = Δ4
    new_evs[5, :, :] = evs[3, :, :]

    return new_evs
end

Δsoln = to_5N_LGE_Δ(Δ_LGE, L=L)
p0 = plot_spatial_profile(Δsoln, L=L, title="Real space")
p1 = plot_in_config_space(Δsoln[5, :, :], L=L, Q=Q, θ=θ, title="On-site")
p2 = plot_in_config_space(Δsoln[1, :, :], L=L, Q=Q, θ=θ, title="-x̂ bond")
p3 = plot_in_config_space(Δsoln[2, :, :], L=L, Q=Q, θ=θ, title="+ŷ bond")
p4 = plot_in_config_space(Δsoln[3, :, :], L=L, Q=Q, θ=θ, title="x̂ bond")
p5 = plot_in_config_space(Δsoln[4, :, :], L=L, Q=Q, θ=θ, title="-ŷ bond")
p = plot(p1, p2, p3, p4, p5, p0, layout=Plots.grid(2, 3,
        widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1000), aspect_ratio=:equal, plot_title=" LGE Δ(J=$J, V0=$V0, V1=$V1, θ=$(θ_to_π(θ)), ϕx=$(θ_to_π(ϕx)), ϕy=$(θ_to_π(ϕy))) for $L × $L lattice at Tc=$(round(Tc,digits=2))")


# Δ = spatial_profile(Δ_LGE, L=L)
# Δ1 = Δ[1, :, :]
# Δ2 = Δ[2, :, :]
# Δ3 = Δ[3, :, :]
# Δ4 = Δ[4, :, :]

# Δ̃1 = circshift(Δ1, (-1, 0)) # mapping between +x̂ and -x̂
# Δ̃2 = circshift(Δ2, (0, 1)) # mapping between +ŷ and -ŷ

# @show maximum(Δ̃1 - Δ3)
# @show maximum(Δ̃2 - Δ4)