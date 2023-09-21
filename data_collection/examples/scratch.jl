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

Js = LinRange(0, 3, 10)
Ts = expspace(-2, 1, 30)

χs_uniform = zeros(length(Js), length(Ts))
χd_uniform = zeros(length(Js), length(Ts))
dχs_uniform = zeros(length(Js), length(Ts))
dχd_uniform = zeros(length(Js), length(Ts))

χs_eigenval = zeros(length(Js), length(Ts))
χd_eigenval = zeros(length(Js), length(Ts))
dχs_eigenval = zeros(length(Js), length(Ts))
dχd_eigenval = zeros(length(Js), length(Ts))

dirpath = joinpath(@__DIR__, "hamiltonians")
mkpath(dirpath)

for (Jᵢ, J) in enumerate(Js)
    loadpath = joinpath(dirpath, "diagonal_$(J)J.h5")
    for (Tᵢ, T) in enumerate(Ts)
        print("$Tᵢ-")

        m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)
        χ_swave, dχ_swave = @time susceptibility_eigenvalue(m, T=T, symmetry="s-wave", calculate_dχdlogT=true, checkpointpath=loadpath)
        χ_dwave, dχ_dwave = @time susceptibility_eigenvalue(m, T=T, symmetry="d-wave", calculate_dχdlogT=true, checkpointpath=loadpath)
        χs_eigenval[Jᵢ, Tᵢ] = χ_swave
        χd_eigenval[Jᵢ, Tᵢ] = χ_dwave
        dχs_eigenval[Jᵢ, Tᵢ] = dχ_swave
        dχd_eigenval[Jᵢ, Tᵢ] = dχ_dwave

        @show dχ_swave, dχ_dwave

        χ, dχdlogT = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true, checkpointpath=loadpath)

        # get the susceptibilities for s-wave and d-wave
        dχ_swave = dχdlogT[1, 1]
        xx, yy, zz = dχdlogT[2, 2], dχdlogT[3, 3], dχdlogT[4, 4]
        xy, yx = dχdlogT[2, 3], dχdlogT[3, 2]
        xz, zx = dχdlogT[2, 4], dχdlogT[4, 2]
        yz, zy = dχdlogT[3, 4], dχdlogT[4, 3]
        dχ_dwave = xx + yy + zz - xy - yx - xz - zx - yz - zy
        @show dχ_swave
        @show dχ_dwave

        dχs_uniform[Jᵢ, Tᵢ] = dχ_swave
        dχd_uniform[Jᵢ, Tᵢ] = dχ_dwave

        # # get the susceptibilities for s-wave and d-wave
        χ_swave = χ[1, 1]
        xx, yy, zz = χ[2, 2], χ[3, 3], χ[4, 4]
        xy, yx = χ[2, 3], χ[3, 2]
        xz, zx = χ[2, 4], χ[4, 2]
        yz, zy = χ[3, 4], χ[4, 3]
        χ_dwave = xx + yy + zz - xy - yx - xz - zx - yz - zy

        χs_uniform[Jᵢ, Tᵢ] = χ_swave
        χd_uniform[Jᵢ, Tᵢ] = χ_dwave
    end
end

# s = copy(dχs_uniform)
# d = copy(dχd_uniform)
# χ_s = copy(χs_uniform)
# χ_d = copy(χd_uniform)

s = copy(dχs_eigenval)
d = copy(dχd_eigenval)
χ_s = copy(χs_eigenval)
χ_d = copy(χd_eigenval)

cmap = cgrad(:matter, length(Ts), categorical=true)
global ps = plot(xlabel="J", ylabel="dχdlogT", title="Eigenvalue s-wave for $L×$L×$L", ylims=(-1.2, 0), size=(800, 800))
global pd = plot(xlabel="J", ylabel="dχdlogT", title="Eigenvalue d-wave for $L×$L×$L", ylims=(-1.2, 0), size=(800, 800))
for (Tᵢ, T) in enumerate(Ts)
    ps = plot!(ps, Js, s[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    pd = plot!(pd, Js, d[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    ps = scatter!(ps, Js, s[:, Tᵢ], c=cmap[Tᵢ], label="T=$T")
    pd = scatter!(pd, Js, d[:, Tᵢ], c=cmap[Tᵢ], label="T=$T")
end

cmap = cgrad(:matter, length(Ts), categorical=true)
global ps = plot(xlabel="J", ylabel="χ", title="Eigenvalue s-wave for $L×$L×$L")
global pd = plot(xlabel="J", ylabel="χ", title="Eigenvalue d-wave for $L×$L×$L")
for (Tᵢ, T) in enumerate(Ts)
    ps = plot!(ps, Js, χ_s[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    pd = plot!(pd, Js, χ_d[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    ps = scatter!(ps, Js, χ_s[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
    pd = scatter!(pd, Js, χ_d[:, Tᵢ], c=cmap[Tᵢ], label=nothing)
end

cmap = cgrad(:matter, length(Js), categorical=true)
global ps = plot(xlabel="T", ylabel="χ", title="Eigenvalue s-wave for $L×$L×$L")
global pd = plot(xlabel="T", ylabel="χ", title="Eigenvalue d-wave for $L×$L×$L")
for (Jᵢ, J) in enumerate(Js)
    ps = plot!(ps, Ts, χ_s[Jᵢ, :], c=cmap[Jᵢ], label=nothing, xaxis=:log10)
    pd = plot!(pd, Ts, χ_d[Jᵢ, :], c=cmap[Jᵢ], label=nothing, xaxis=:log10)
    ps = scatter!(ps, Ts, χ_s[Jᵢ, :], c=cmap[Jᵢ], label="J=$J", xaxis=:log10)
    pd = scatter!(pd, Ts, χ_d[Jᵢ, :], c=cmap[Jᵢ], label="J=$J", xaxis=:log10)
end

cmap = cgrad(:matter, length(Js), categorical=true)
global ps = plot(xlabel="T", ylabel="dχdlogT", title="Eigenvalue s-wave for $L×$L×$L")
global pd = plot(xlabel="T", ylabel="dχdlogT", title="Eigenvalue d-wave for $L×$L×$L")
for (Jᵢ, J) in enumerate(Js)
    ps = plot!(ps, Ts, s[Jᵢ, :], c=cmap[Jᵢ], label=nothing, xaxis=:log10)
    pd = plot!(pd, Ts, d[Jᵢ, :], c=cmap[Jᵢ], label=nothing, xaxis=:log10)
    ps = scatter!(ps, Ts, s[Jᵢ, :], c=cmap[Jᵢ], label="J=$J", xaxis=:log10)
    pd = scatter!(pd, Ts, d[Jᵢ, :], c=cmap[Jᵢ], label=nothing, xaxis=:log10)
end