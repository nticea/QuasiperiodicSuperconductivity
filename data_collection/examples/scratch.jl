## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

# Parameters 
t = 1
L = 7
Q = (√5 - 1) / 2
μ = 1
θ = π / 7
V0 = 0
V1 = 0
ϕx = 0#rand() * 2π
ϕy = 0#rand() * 2π
ϕz = 0#rand() * 2π
periodic = false
disorder = false
ndims = 3

Λ = nothing

Ts = expspace(-3, 1, 20)
Js = [0, 0.2, 4]
χs_swave, dχs_swave = zeros(length(Ts), length(Js)), zeros(length(Ts), length(Js))
χs_dwave, dχs_dwave = zeros(length(Ts), length(Js)), zeros(length(Ts), length(Js))

## CALCULATION ## 
for (Jᵢ, J) in enumerate(Js)
    for (Tᵢ, T) in enumerate(Ts)
        m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        χ, dχdlogT = @time uniform_susceptibility(m, T=T, calculate_dχdlogT=true, Λ=Λ)
        χswave, χdwave = uniform_susceptibility_components(χ, ndims=ndims)
        dχswave, dχdwave = uniform_susceptibility_components(dχdlogT, ndims=ndims)
        χs_swave[Tᵢ, Jᵢ] = χswave
        χs_dwave[Tᵢ, Jᵢ] = χdwave
        dχs_swave[Tᵢ, Jᵢ] = dχswave
        dχs_dwave[Tᵢ, Jᵢ] = dχdwave
    end
end

p = plot(xlabel="T", ylabel="dχ/dlogT")
cmap = cgrad(:matter, length(Js), categorical=true)
for (Jᵢ, J) in enumerate(Js)
    plot!(p, Ts, dχs_dwave[:, Jᵢ], xaxis=:log10, label=nothing, c=cmap[Jᵢ])
    scatter!(p, Ts, dχs_dwave[:, Jᵢ], xaxis=:log10, label="J=$J", c=cmap[Jᵢ])
end


