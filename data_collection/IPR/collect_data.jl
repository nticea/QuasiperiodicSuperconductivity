## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")
include("../../src/utilities.jl")

## MODEL PARAMETERS ##
args = parse.(Float64, ARGS)
@assert length(args) == 14
L, t, Q, μ, θ, ϕx, ϕy, ϕz, V0, V1, J, periodic, ndims, disorder = [args[n] for n in 1:length(args)]
L = Int(L)
ndims = Int(ndims)
periodic = Bool(periodic)
disorder = Bool(disorder)
m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)

if disorder
    pot = "disorder"
else
    pot = "QP"
end
# datapath
scratchbase = joinpath("/scratch/users/nticea", "QuasiperiodicSuperconductivity", "IPR")
dirname = "$(ndims)D_data_$(pot)"
datapath = joinpath(scratchbase, dirname)
mkpath(datapath)

df = DataFrame(L=[], μ=[], J=[], ϕx=[], ϕy=[], ϕz=[], ipr_real=[], ipr_k=[], E=[], pot=[])

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
H = DiagonalizedHamiltonian(m)
E = H.E
# calculate the real IPR 
ipr_real = IPR_real(H)
# calculate the momentum IPR 
ipr_k = IPR_momentum(H)

timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
savepath = joinpath(datapath, "$(L)L_$(J)J" * timestamp * ".csv")

dfi = DataFrame(L=[L], μ=[μ], J=[J], ϕx=[ϕx], ϕy=[ϕy], ϕz=[ϕz], ipr_real=[ipr_real], ipr_k=[ipr_k], E=[E], pot=[pot])
append!(df, dfi)

if isfile(savepath)
    CSV.write(savepath, df, append=true)
else
    CSV.write(savepath, df, append=false)
end

