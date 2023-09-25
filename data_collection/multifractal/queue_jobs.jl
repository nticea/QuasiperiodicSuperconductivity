## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("../../src/submit_job.jl")

## PARAMETERS ## 
Js = LinRange(0, 5, 30)
ℓ = 5
E₀ = 0
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = -1
V1 = -1
ϕx = 0
ϕy = 0
ϕz = 0
periodic = 0
ndims = 3

filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "multifractal_scaling"

for J in Js
    # L = 10
    ps = ModelParams(L=10, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=100, kwargs="$ℓ $E₀", time="7:00")

    # L = 15
    ps = ModelParams(L=15, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=200, kwargs="$ℓ $E₀", time="1:00:00")

    # L = 20
    ps = ModelParams(L=20, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=300, kwargs="$ℓ $E₀", time="3:00:00")

    # L = 25
    ps = ModelParams(L=25, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=350, kwargs="$ℓ $E₀", time="5:00:00")
end
