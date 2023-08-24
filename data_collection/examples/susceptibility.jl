
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/meanfield.jl")

## PARAMETERS ## 
L = 23
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
ϕx = 0
ϕy = 0
ϕz = 0
periodic = true
ndims = 3

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
Ts = expspace(-2, 1.3, 20)
χswave = []
χdwave = []

for T in Ts
    @show T
    χ = uniform_susceptibility(m, T=T, symmetry="d-wave")
    xx, yy = χ[2, 2], χ[3, 3]
    xy, yx = χ[2, 3], χ[3, 2]
    if ndims == 2
        d = xx + yy - xy - yx
    elseif ndims == 3
        zz = χ[4, 4]
        xz, zx = χ[2, 4], χ[4, 2]
        yz, zy = χ[3, 4], χ[4, 3]
        d = xx + yy + zz - xy - yx - xz - zx - yz - zy
    else
        println("sorry")
    end

    s = χ[1, 1]

    push!(χswave, s)
    push!(χdwave, d)
end

plot(Ts, χdwave, xaxis=:log10)
plot(Ts, χswave, xaxis=:log10)