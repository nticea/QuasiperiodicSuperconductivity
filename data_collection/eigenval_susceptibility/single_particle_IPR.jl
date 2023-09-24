## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")
include("../../src/results.jl")

## PARAMETERS ## 
L = 11
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 1.5
V1 = -1
ϕxs = LinRange(0, π + 0.01, 3)
ϕys = LinRange(0, π + 0.1, 3)
ϕzs = LinRange(0, π - 0.01, 3)
periodic = true
ndims = 3
nbins = 41

Js = collect(0:0.25:3)
dos = zeros(length(Js), nbins)

df = DataFrame(J=[], ϕx=[], ϕy=[], ϕz=[], ipr_real=[], ipr_k=[])

for (j, J) in enumerate(Js)
    print(j, "-")
    for (ix, ϕx) in enumerate(ϕxs)
        for (iy, ϕy) in enumerate(ϕys)
            for (iz, ϕz) in enumerate(ϕzs)
                # initialize model 
                m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
                H = DiagonalizedHamiltonian(m)

                # calculate the real IPR 
                ipr = IPR_real(H)
                ipr_bin = bin_results(ipr, nbins=nbins)
                # take the middle bin 
                ipr_real = mean(ipr)#ipr_bin[floor(Int, nbins / 2)]

                # calculate the real IPR 
                ipr = IPR_momentum(H)
                ipr_bin = bin_results(ipr, nbins=nbins)
                # take the middle bin
                ipr_k = mean(ipr_bin)#ipr_bin[floor(Int, nbins / 2)]
                @show mean(ipr)

                dfi = DataFrame(J=[J], ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ipr_real=[ipr_real], ipr_k=[ipr_k])
                append!(df, dfi)
            end
        end
    end
end

# use split-apply-combine
gdf = groupby(df, [:J])
dfmean = combine(gdf, [:ipr_real => mean, :ipr_k => mean])

plot(dfmean.J, dfmean.ipr_real_mean, color="red")
scatter!(dfmean.J, dfmean.ipr_real_mean, color="red")
plot!(dfmean.J, dfmean.ipr_k_mean, color="blue")
scatter!(dfmean.J, dfmean.ipr_k_mean, color="blue")