## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")
include("../../src/results.jl")
include("../../src/utilities.jl")

## PARAMETERS ## 
Ls = [18, 27, 30]
Js = LinRange(1, 3, 20)
ℓ = 3
E₀ = 1
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = -1
V1 = -1
ϕx = rand(0:2π, 5)
ϕy = rand(0:2π, 5)
ϕz = 0
periodic = true
ndims = 2
nbins = 31

for (l, L) in enumerate(Ls)
    print("$l-")
    for (j, J) in enumerate(Js)
        m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
        α₀ = multifractal_mean(m; E₀=E₀)

        # save it out 
        save_results(m, ℓ=ℓ, E₀=E₀, α₀=α₀)
    end
end

# p = plot()
# cmap = ["red", "blue", "green", "orange", "black"]
# for (l, L) in enumerate(Ls)
#     p = plot!(Js, αs[l, :], c=cmap[l])
#     p = scatter!(p, Js, αs[l, :], c=cmap[l], label="L=$L")
# end
