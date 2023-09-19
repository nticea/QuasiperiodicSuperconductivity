## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

# Parameters 
L = 7
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

Js = LinRange(0, 3, 12)
Ts = expspace(-2, 1, 30)

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)

s = zeros(length(Js), length(Ts))
d = zeros(length(Js), length(Ts))
for (Tᵢ, T) in enumerate(Ts)
    for (Jᵢ, J) in enumerate(Js)
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

        s[Jᵢ, Tᵢ] = dχ_swave
        d[Jᵢ, Tᵢ] = dχ_dwave
    end
end

# ps = plot(Ts, s, label=nothing, xlabel="T", ylabel="dχdlogT", xaxis=:log10)
# pd = plot(Ts, d, label=nothing, xlabel="T", ylabel="dχdlogT", xaxis=:log10)
# ps = scatter!(ps, Ts, s, label="s-wave", xlabel="T", ylabel="dχdlogT", xaxis=:log10)
# pd = scatter!(pd, Ts, d, label="d-wave", xlabel="T", ylabel="dχdlogT", xaxis=:log10)
