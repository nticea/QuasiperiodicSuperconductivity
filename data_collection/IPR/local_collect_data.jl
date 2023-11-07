## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")
include("../../src/utilities.jl")

# Parameters 
L = 11
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
periodic = true
ndims = 3

dirpath = joinpath(@__DIR__, "data")
mkpath(dirpath)
savepath = joinpath(dirpath, "IPR_data_$(L)L.csv")
df = DataFrame(L=[], μ=[], J=[], ϕx=[], ϕy=[], ϕz=[], ipr_real=[], ipr_k=[], E=[], pot=[])

Js = collect(0:0.25:10)
nrep = 1

# QUASIPERIODIC
for n in 1:nrep
    ϕx, ϕy, ϕz = 2π * rand(), 2π * rand(), 2π * rand()
    for J in Js
        print("n=$n,J=$J-")
        m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
        H = DiagonalizedHamiltonian(m)
        E = H.E
        # calculate the real IPR 
        ipr_real = IPR_real(H)
        # calculate the momentum IPR 
        ipr_k = IPR_momentum(H)

        dfi = DataFrame(L=[L], μ=[μ], J=[J], ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ipr_real=[ipr_real], ipr_k=[ipr_k], E=[E], pot=["QP"])
        append!(df, dfi)

        if isfile(savepath)
            CSV.write(savepath, df, append=true)
        else
            CSV.write(savepath, df, append=false)
        end
    end
end

# DISORDER 
for n in 1:nrep
    ϕx, ϕy, ϕz = 2π * rand(), 2π * rand(), 2π * rand()
    for J in Js
        print("n=$n,J=$J-")
        m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
        H = DiagonalizedHamiltonian(m)
        E = H.E
        # calculate the real IPR 
        ipr_real = IPR_real(H)
        # calculate the momentum IPR 
        ipr_k = IPR_momentum(H)

        dfi = DataFrame(L=[L], μ=[μ], J=[J], ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ipr_real=[ipr_real], ipr_k=[ipr_k], E=[E], pot=["disorder"])
        append!(df, dfi)

        if isfile(savepath)
            CSV.write(savepath, df, append=true)
        else
            CSV.write(savepath, df, append=false)
        end
    end
end
