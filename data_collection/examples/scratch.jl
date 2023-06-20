
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
ϕx = 0.239
ϕy = 0.547
V0 = 1
V1 = -1.5
J = 1
periodic = true
niter = 500
tol = 1e-12

T = 0.1
λ_LGE, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry="d-wave")

# p_BdG = plot_spatial_profile(Δ_BdG ./ (γ * λ_LGE); L, title="BdG Δ after one iteration")
# p_LGE = plot_spatial_profile(Δ_LGE; L, title="LGE Δ (initial)")
# p = plot(p_LGE, p_BdG, layout=Plots.grid(1, 2, widths=[1 / 2, 1 / 2]), size=(1500, 1000), aspect_ratio=:equal)

#@assert 1 == 0

γs = LinRange(1e-9, 1, 50)
# max_ΔBdG = []
# max_ΔLGE = []
diff = []
for (n, γ) in enumerate(γs)
    print(n, "-")
    Δ_BdG, Δmax = Δ_dwave_debug(T, γ, λ_LGE; L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, niter=niter, tol=tol, Δ_init=Δ_LGE)
    push!(diff, norm(Δ_BdG - λ_LGE * γ * Δ_LGE))
    # push!(max_ΔBdG, maximum(abs.(Δ_BdG)) / γ)
    # push!(max_ΔLGE, maximum(abs.(Δ_LGE)) * λ_LGE)
    #push!(maxs, Δmax)
end

# plot(γs, maxs, label=nothing, c="blue")
# scatter!(γs, maxs, label=nothing, c="blue")
# hline!([λ_LGE], label="λ LGE")

plot(γs, diff, c="blue")
#plot!(γs, max_ΔLGE, label="LGE × λ_LGE", c="red")
xaxis!("γ")
#yaxis!("||Δ' - γ×λ×ΔLGE||/γ")
#title!("Checking BdG")
