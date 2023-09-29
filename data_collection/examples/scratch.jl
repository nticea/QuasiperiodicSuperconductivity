## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

# Parameters 
t = 1
J = 0.25
Q = (√5 - 1) / 2
μ = 1
θ = π / 7
V0 = 0
V1 = 0
ϕx = 0
ϕy = 0
ϕz = 0
periodic = false
ndims = 3

Λ = 0.3

Ts = expspace(-2, 1, 20)
Ls = [6, 7]
χs_swave, dχs_swave = zeros(length(Ts), length(Ls)), zeros(length(Ts), length(Ls))
χs_dwave, dχs_dwave = zeros(length(Ts), length(Ls)), zeros(length(Ts), length(Ls))

## CALCULATION ## 
for (Lᵢ, L) in enumerate(Ls)
    for (Tᵢ, T) in enumerate(Ts)
        m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
        χ, dχdlogT = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true)
        χswave, χdwave = uniform_susceptibility_components(χ, ndims=ndims)
        dχswave, dχdwave = uniform_susceptibility_components(dχdlogT, ndims=ndims)
        χs_swave[Tᵢ, Lᵢ] = χswave
        χs_dwave[Tᵢ, Lᵢ] = χdwave
        dχs_swave[Tᵢ, Lᵢ] = dχswave
        dχs_dwave[Tᵢ, Lᵢ] = dχdwave
        @show J, T, χ, dχdlogT
    end
end

p = plot(xlabel="T", ylabel="dχ/dlogT")
cmap = cgrad(:matter, length(Ls), categorical=true)
for (Lᵢ, L) in enumerate(Ls)
    p = plot!(p, Ts, dχs_dwave[:, Lᵢ], xaxis=:log10, label=nothing, c=cmap[Lᵢ])
    p = scatter!(p, Ts, dχs_dwave[:, Lᵢ], xaxis=:log10, label="L=$L", c=cmap[Lᵢ])
end


