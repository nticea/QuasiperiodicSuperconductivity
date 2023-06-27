
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
include("../../src/BdG.jl")

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
J = 0.5
periodic = true
pairing_symmetry = "d-wave"
niter = 500
tol = 1e-12

T = 0
λ_LGE, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")

@error "Check that the energy is monotonically decreasing as we iterate using expectation value of H0 + 1/2 H_int"
@error "Check whether BdG is flowing away from half-filling by computing the density"

γs = LinRange(1e-9, 1, 50)
diff = []
for (n, γ) in enumerate(γs)
    print(n, "-")
    Δ_BdG, hist, Es = Δ_dwave_debug(T, γ, λ_LGE; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, Δ_init=Δ_LGE)
    push!(diff, norm(Δ_BdG - λ_LGE * γ * Δ_LGE))
end

plot(γs, diff, c="blue")
xaxis!("γ")