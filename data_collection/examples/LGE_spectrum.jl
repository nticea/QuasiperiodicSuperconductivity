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
L = 11 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0.375
ϕy = 0.29
V0 = 1.5
V1 = -1
periodic = true
pairing_symmetry = "d-wave"

T = 0.04
Js = collect(0:0.1:2)
if pairing_symmetry == "d-wave"
    N = 5 * L * L
elseif pairing_symmetry == "s-wave"
    N = L * L
else
    @error "Invalid pairing symmetry"
end
spectrum = zeros(N, length(Js))

# generate data 
for (i, J) in enumerate(Js)
    vals = LGE_spectrum(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)
    spectrum[:, i] = vals
end

# plot everything
cmap = reverse(cgrad(:matter, N, categorical=true))
p1 = plot()
for i in 1:N
    vals = spectrum[i, :]
    plot!(p1, Js, vals, c=cmap[i], label=nothing)
end
xlabel!(p1, "J")
title!("LGE spectrum for $pairing_symmetry pairing")
plot(p1)

