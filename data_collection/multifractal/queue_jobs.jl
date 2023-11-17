## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("../../src/submit_job.jl")

## PARAMETERS ## 
J = 0
λ = 1 / 5
E₀ = 0.75
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

nrep = 30

filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "multifractal_scaling"

for _ in 1:nrep
    # make random potential offsets 
    ϕx = 2π * rand()
    ϕy = 2π * rand()
    ϕz = 2π * rand()

    # L = 10
    ps = ModelParams(L=10, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=32, kwargs="$λ $E₀", time="5:00")

    # L = 15
    ps = ModelParams(L=15, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=64, kwargs="$λ $E₀", time="15:00")

    # L = 20
    ps = ModelParams(L=20, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=90, kwargs="$λ $E₀", time="30:00")

    # L = 25
    ps = ModelParams(L=25, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, kwargs="$λ $E₀", time="2:00:00")
end
