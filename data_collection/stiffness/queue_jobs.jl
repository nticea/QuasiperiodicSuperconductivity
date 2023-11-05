## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("../../src/submit_job.jl")

t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0.75
θ = π / 7
periodic = 1
ndims = 3
disorder = false

Js = collect(0:0.25:3)
filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "stiffness"

nrep = 40

for _ in 1:nrep

    # make random potential offsets 
    ϕx = 2π * rand()
    ϕy = 2π * rand()
    ϕz = 2π * rand()

    for J in Js

        # d-wave 
        V0 = 2
        V1 = -1.5
        # ps = ModelParams(L=13, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        # submit_job(ps, filepath, @__DIR__, job_prefix, mem=256, time="6:00:00")
        ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        submit_job(ps, filepath, @__DIR__, job_prefix, mem=256, time="4:00:00")
        # ps = ModelParams(L=7, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        # submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="1:00:00")

        # s-wave 
        V0 = -3
        V1 = 0
        # ps = ModelParams(L=13, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        # submit_job(ps, filepath, @__DIR__, job_prefix, mem=256, time="6:00:00")
        ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        submit_job(ps, filepath, @__DIR__, job_prefix, mem=256, time="4:00:00")
        # ps = ModelParams(L=7, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
        # submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="1:00:00")
    end
end
