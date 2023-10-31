## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames
using Dates

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")
include("utilities.jl")

## MODEL PARAMETERS ##
args = parse.(Float64, ARGS)
@assert length(args) == 14
L, t, Q, μ, θ, ϕx, ϕy, ϕz, V0, V1, J, periodic, ndims, disorder = [args[n] for n in 1:length(args)]
L = Int(L)
ndims = Int(ndims)
periodic = Bool(periodic)
disorder = Bool(disorder)
m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)

@show V0, V1, J

## SIMULATION PARAMETERS ## 
niter = 500
BdG_tol = 1e-15
LGE_tol = 1e-2

## SAVING ## 
mkpath(joinpath(@__DIR__, "data"))
dfs = load_dfs()
timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
savepath_BdG = joinpath(@__DIR__, "data", "$(L)L_ΦQ_BdG_" * timestamp * ".csv")
savepath_LGE = joinpath(@__DIR__, "data", "$(L)L_ΦQ_LGE_" * timestamp * ".csv")
df_BdG = DataFrame(L=[], t=[], μ=[], J=[], Q=[], θ=[],
    ϕx=[], ϕy=[], ϕz=[], V0=[], V1=[], ndims=[], periodic=[],
    T=[], λ=[], Δ=[], K=[], Π=[])
df_LGE = DataFrame(L=[], t=[], μ=[], J=[], Q=[], θ=[],
    ϕx=[], ϕy=[], ϕz=[], V0=[], V1=[], ndims=[], periodic=[],
    T=[], λ=[], Δ=[], K=[], Π=[])

if !already_computed(m, dfs, T=0)

    ## Tc using LGE ##
    println("Finding Tc using LGE")
    Tc, λ, Δ_LGE = LGE_find_Tc(m, npts=5, tol=LGE_tol)
    update_results!(m, df_LGE; T=Tc, λ=λ, Δ=Δ_LGE)
    CSV.write(savepath_LGE, df_LGE)
    flush(stdout)

    ## Superfluid stiffness calculation ##
    if isfinite(Tc)

        # Everything is at 0T
        T = 0

        # Get the initial LGE guess 
        println("Finding LGE sol'n at T=0")
        λ, Δ_LGE = @time pairfield_correlation(m, T=T)
        update_results!(m, df_LGE; T=T, λ=λ, Δ=Δ_LGE)
        CSV.write(savepath_LGE, df_LGE)
        flush(stdout)

        # Superfluid stiffness
        println("Computing superfluid stiffness at T=0")
        K, Π, Δ_BdG = @time superfluid_stiffness_finiteT(m, T=T, tol=BdG_tol, niter=niter, Δ_init=Δ_LGE)
        update_results!(m, df_BdG; T=T, λ=λ, Δ=Δ_BdG, K=K, Π=Π)
        CSV.write(savepath_BdG, df_BdG)
        flush(stdout)
    end
end