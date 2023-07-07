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
L = 7 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
ϕx = 0.1234
ϕy = 0.5678
V0 = 1
V1 = -1.5
periodic = true
niter = 500
tol = 1e-15
T = 0

J = 0.5
λ, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)
K, Π = superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, Δ_init=Δ_LGE)

# @error "Please check factors of 2 in the kinetic, Aq, and Dq terms!"

# Js = collect(0:0.1:2)
# Ks = zeros(length(Js), 4)
# Πs = zeros(length(Js), 4)

# for (j, J) in enumerate(Js)
#     @show j
#     λ, Δ_LGE = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)
#     K, Π = superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, Δ_init=Δ_LGE)

#     Ks[j, :] = K
#     Πs[j, :] = Π
# end

# Js = collect(0:0.1:2)
# Tcs_dwave = zeros(length(Js))
# for (j, J) in enumerate(Js)
#     Tcs_dwave[j] = LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, tol=1e-3)
# end

# Ds = (-Ks + Πs) / π

# p = plot()
# dirs = ["x̂", "ŷ", "-x̂", "-ŷ"]
# cmap = cgrad(:matter, 4, categorical=true)

# for i in 1:4
#     plot!(p, Js, Ds[:, i], label=nothing, c=cmap[i])
#     scatter!(p, Js, Ds[:, i], label=dirs[i], c=cmap[i])
# end
# title!(p, "Superfluid stiffness for \n Δ(V0=$V0, V1=$V1, θ=$(θ_to_π(θ)),ϕx=$ϕx,ϕy=$ϕy) \n $L × $L lattice at T=$(round(T,digits=2))")
# xlabel!("J")
# ylabel!("Tc")

# plot!(Js, Tcs_dwave, label=nothing, c="blue")
# scatter!(Js, Tcs_dwave, label="LGE soln", c="blue")


