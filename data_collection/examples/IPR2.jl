## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")
include("../../src/utilities.jl")

# Parameters 
L = 13
λ = 1 / 13
E₀ = 0.75
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
periodic = true
disorder = false
ndims = 3

if periodic
    BC = "PBC"
else
    BC = "OBC"
end
if disorder
    pot = "disorder"
else
    pot = "QP"
end

savepath = "IPR_data_$(BC)_$(pot).csv"

df = DataFrame(L=[], J=[], ϕx=[], ϕy=[], ϕz=[], ipr_real=[], ipr_k=[], α₀=[], BC=[], pot=[])
# fname = "/Users/nicole/Dropbox/Grad/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/IPR_data.csv"
# df = DataFrame(CSV.File(fname))

Js = collect(0:0.25:10)
nrep = 20
for n in 1:nrep
    ϕx, ϕy, ϕz = 2π * rand(), 2π * rand(), 2π * rand()
    for J in Js
        print("n=$n,J=$J-")
        m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        α₀, ipr_real, ipr_k = compute_scaling_properties(m; λ=λ, E₀=E₀)

        dfi = DataFrame(L=[L], J=[J], ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ipr_real=[ipr_real], ipr_k=[ipr_k], α₀=[real(α₀)], BC=[BC], pot=[pot])
        append!(df, dfi)

        if isfile(savepath)
            CSV.write(savepath, df, append=true)
        else
            CSV.write(savepath, df, append=false)
        end
    end
end

Js = collect(0:0.25:6)
nrep = 20
for n in 1:nrep
    ϕx, ϕy, ϕz = 2π * rand(), 2π * rand(), 2π * rand()
    for J in Js
        print("n=$n,J=$J-")
        m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        α₀, ipr_real, ipr_k = compute_scaling_properties(m; λ=λ, E₀=E₀)

        dfi = DataFrame(L=[L], J=[J], ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ipr_real=[ipr_real], ipr_k=[ipr_k], α₀=[real(α₀)], BC=[BC], pot=[pot])
        append!(df, dfi)

        if isfile(savepath)
            CSV.write(savepath, df, append=true)
        else
            CSV.write(savepath, df, append=false)
        end
    end
end

# use split-apply-combine
df_plot = df[(df.BC.==BC).&(df.pot.==pot), :]
gdf = groupby(df_plot, [:L, :J])
dfmean = combine(gdf, [:ipr_real => mean, :ipr_k => mean])
dfmean = dfmean[(dfmean.L.==L), :]

plot(dfmean.J, dfmean.ipr_real_mean, color="red", label=nothing)
scatter!(dfmean.J, dfmean.ipr_real_mean, color="red", label="real space")
plot!(dfmean.J, dfmean.ipr_k_mean, color="blue", label=nothing)
scatter!(dfmean.J, dfmean.ipr_k_mean, color="blue", label="momentum space")
title!("IPR ($pot $L×$L×$L at E₀=$E₀, μ=$μ with $BC)")
xlabel!("J")
ylabel!("IPR")
