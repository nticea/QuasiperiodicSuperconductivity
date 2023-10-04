## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")
include("../../src/results.jl")

## PARAMETERS ## 
Ls = [5, 10]
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 1.5
V1 = -1
ϕxs = LinRange(0, π + 0.01, 3)
ϕys = LinRange(0, π + 0.1, 3)
ϕzs = LinRange(0, π - 0.01, 3)
periodic = false
ndims = 3

Js = collect(0:0.25:6)
df = DataFrame(L=[], J=[], ϕx=[], ϕy=[], ϕz=[], ipr_real=[], ipr_k=[])

# m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=0, ϕy=0, ϕz=0, V0=V0, V1=V1, J=0, periodic=periodic, ndims=ndims)
# H = DiagonalizedHamiltonian(m)
# Uq = fourier_transform_eigenstates(H)

# # IPR = ∑ₖ |ψ(k)|⁴ 
# u = Uq[:, 666]
# u2 = u .* conj.(u)
# u4 = u2 .^ 2

# @assert 1 == 0

for (l, L) in enumerate(Ls)
    for (j, J) in enumerate(Js)
        print("J=$J,L=$L-")
        for (ix, ϕx) in enumerate(ϕxs)
            for (iy, ϕy) in enumerate(ϕys)
                for (iz, ϕz) in enumerate(ϕzs)
                    # initialize model 
                    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
                    H = DiagonalizedHamiltonian(m)

                    # calculate the real IPR 
                    ipr = IPR_real(H)
                    mididx = floor(Int, length(ipr) / 2)
                    ipr_real = ipr[mididx]
                    @show ipr_real

                    # calculate the real IPR 
                    ipr = IPR_momentum(H)
                    ipr_k = ipr[mididx]
                    @show ipr_k

                    dfi = DataFrame(L=[L], J=[J], ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ipr_real=[ipr_real], ipr_k=[ipr_k])
                    append!(df, dfi)
                end
            end
        end
    end
end

# use split-apply-combine
L=10
gdf = groupby(df, [:L, :J])
dfmean = combine(gdf, [:ipr_real => mean, :ipr_k => mean])
dfmean = dfmean[(dfmean.L.==L), :]

plot(dfmean.J, dfmean.ipr_real_mean, color="red", label=nothing)
scatter!(dfmean.J, dfmean.ipr_real_mean, color="red", label="real space")
plot!(dfmean.J, dfmean.ipr_k_mean, color="blue", label=nothing)
scatter!(dfmean.J, dfmean.ipr_k_mean, color="blue", label="momentum space")
title!("IPR")
xlabel!("J")
ylabel!("IPR")