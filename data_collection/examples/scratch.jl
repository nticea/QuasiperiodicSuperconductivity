## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

# Parameters 
L = 17
t = 1
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

J = 0.3
T = 1e-2

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)
χ, dχdlogT = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true)

# get the susceptibilities for s-wave and d-wave
dχ_swave = dχdlogT[1, 1]
xx, yy, zz = dχdlogT[2, 2], dχdlogT[3, 3], dχdlogT[4, 4]
xy, yx = dχdlogT[2, 3], dχdlogT[3, 2]
xz, zx = dχdlogT[2, 4], dχdlogT[4, 2]
yz, zy = dχdlogT[3, 4], dχdlogT[4, 3]
dχ_dwave = xx + yy + zz - xy - yx - xz - zx - yz - zy

@show dχ_swave
@show dχ_dwave