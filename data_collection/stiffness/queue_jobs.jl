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

Js = collect(0:0.25:6)
filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "stiffness"

# make random potential offsets hi
ϕx = 0
ϕy = 0
ϕz = 0

for J in Js
    # QUASIPERIODIC
    # d-wave 
    V0 = 1
    V1 = -3
    ps = ModelParams(L=13, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=800, time="12:00:00")
    ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="8:00:00")
    ps = ModelParams(L=7, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="1:00:00")

    # s-wave 
    V0 = -3
    V1 = 0
    ps = ModelParams(L=13, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=800, time="12:00:00")
    ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="8:00:00")
    ps = ModelParams(L=7, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="1:00:00")

    # DISORDER 
    # d-wave 
    V0 = 1
    V1 = -3
    ps = ModelParams(L=13, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=800, time="12:00:00")
    ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="8:00:00")
    ps = ModelParams(L=7, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="1:00:00")

    # s-wave 
    V0 = -3
    V1 = 0
    ps = ModelParams(L=13, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=800, time="12:00:00")
    ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="8:00:00")
    ps = ModelParams(L=7, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="1:00:00")
end

