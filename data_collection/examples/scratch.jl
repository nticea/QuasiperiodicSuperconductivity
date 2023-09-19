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

s = zeros(length(Js), length(Ts))
d = zeros(length(Js), length(Ts))
χ_s = zeros(length(Js), length(Ts))
χ_d = zeros(length(Js), length(Ts))
for (Tᵢ, T) in enumerate(Ts)
    for (Jᵢ, J) in enumerate(Js)
        println("$J-")

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

        s[Jᵢ, Tᵢ] = dχ_swave
        d[Jᵢ, Tᵢ] = dχ_dwave

        # get the susceptibilities for s-wave and d-wave
        χ_swave = χ[1, 1]
        xx, yy, zz = χ[2, 2], χ[3, 3], χ[4, 4]
        xy, yx = χ[2, 3], χ[3, 2]
        xz, zx = χ[2, 4], χ[4, 2]
        yz, zy = χ[3, 4], χ[4, 3]
        χ_dwave = xx + yy + zz - xy - yx - xz - zx - yz - zy

        χ_s[Jᵢ, Tᵢ] = χ_swave
        χ_d[Jᵢ, Tᵢ] = χ_dwave
    end
end

# ps = plot(Ts, s, label=nothing, xlabel="T", ylabel="dχdlogT", xaxis=:log10)
# pd = plot(Ts, d, label=nothing, xlabel="T", ylabel="dχdlogT", xaxis=:log10)
# ps = scatter!(ps, Ts, s, label="s-wave", xlabel="T", ylabel="dχdlogT", xaxis=:log10)
# pd = scatter!(pd, Ts, d, label="d-wave", xlabel="T", ylabel="dχdlogT", xaxis=:log10)

cmap = cgrad(:matter, length(Ts), categorical=true)
ps = plot(xlabel="J", ylabel="dχdlogT", title="s-wave")
pd = plot(xlabel="J", ylabel="dχdlogT", title="s-wave")
for (Tᵢ, T) in enumerate(Ts)
    ps = plot!(ps, Js, s[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    pd = plot!(pd, Js, d[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    ps = scatter!(ps, Js, s[:, Tᵢ], c=cmap[Tᵢ], label="T=$T")
    pd = scatter!(pd, Js, d[:, Tᵢ], c=cmap[Tᵢ], label="T=$T")
end

cmap = cgrad(:matter, length(Ts), categorical=true)
ps = plot(xlabel="J", ylabel="χ", title="s-wave")
pd = plot(xlabel="J", ylabel="χ", title="s-wave")
for (Tᵢ, T) in enumerate(Ts)
    ps = plot!(ps, Js, χ_s[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    pd = plot!(pd, Js, χ_d[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    ps = scatter!(ps, Js, χ_s[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    pd = scatter!(pd, Js, χ_d[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
end