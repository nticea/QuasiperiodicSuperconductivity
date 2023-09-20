## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("../../src/submit_job.jl")

## PARAMETERS ## 
Js = LinRange(0, 5, 30)
ℓ = 3
E₀ = 1
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = -1
V1 = -1
ϕxs = LinRange(0, 2π, 3)
ϕys = LinRange(0, 2π, 3)
ϕz = 0
periodic = 1
ndims = 3

filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "multifractal_scaling"

for ϕx in ϕxs
    for ϕy in ϕys
        for J in Js
            # L = 18
            ps = ModelParams(L=18, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
            submit_job(ps, filepath, @__DIR__, job_prefix, mem=50, kwargs="$ℓ $E₀", time="5:00")

            # L = 27
            ps = ModelParams(L=27, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
            submit_job(ps, filepath, @__DIR__, job_prefix, mem=250, kwargs="$ℓ $E₀", time="4:00:00")

            # L = 30
            ps = ModelParams(L=30, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
            submit_job(ps, filepath, @__DIR__, job_prefix, mem=300, kwargs="$ℓ $E₀", time="6:00:00")
        end
    end
end