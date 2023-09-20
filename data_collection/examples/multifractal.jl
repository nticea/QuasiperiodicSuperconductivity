
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")

## PARAMETERS ## 
L = 11
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 1.5
V1 = -1
ϕx = 0
ϕy = 0
ϕz = 0
periodic = true
ndims = 3
E₀ = 1
ℓ = 1

Js = LinRange(0, 5, 30)
αs = []
for (Jᵢ, J) in enumerate(Js)
    print("$Jᵢ-")
    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)
    α₀ = multifractal_mean(m; E₀=E₀, ℓ=ℓ)
    push!(αs, α₀)
    @show α₀
end