
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots, Statistics
include("../../src/model.jl")

## PARAMETERS ## 
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 1.5
V1 = -1
ϕxs = LinRange(0, π, 3)
ϕys = LinRange(0, π, 3)
ϕz = 0
periodic = false
ndims = 3
E₀ = 1
ℓ = 1

Js = LinRange(2, 4, 20)
Ls = [5, 10, 15]
αs = zeros(length(Ls), length(Js), length(ϕxs), length(ϕys))
for (xᵢ, ϕx) in enumerate(ϕxs)
    for (yᵢ, ϕy) in enumerate(ϕys)
        for (Lᵢ, L) in enumerate(Ls)
            for (Jᵢ, J) in enumerate(Js)
                print("$Jᵢ-")
                m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)
                α₀ = multifractal_mean(m; E₀=E₀, ℓ=ℓ)
                αs[Lᵢ, Jᵢ, xᵢ, yᵢ] = α₀
                @show α₀
            end
        end
    end
end

αs_mean = mean(αs, dims=[3, 4])

cmap = cgrad(:matter, length(Ls), categorical=true)
p = plot(xlabel="J", ylabel="α₀")
for (Lᵢ, L) in enumerate(Ls)
    p = plot!(p, Js, αs_mean[Lᵢ, :], label=nothing, c=cmap[Lᵢ])
    p = scatter!(Js, αs_mean[Lᵢ, :], label="L=$L", c=cmap[Lᵢ])
end