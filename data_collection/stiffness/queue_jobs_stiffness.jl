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
filepath = joinpath(@__DIR__, "collect_data_stiffness.jl")
job_prefix = "stiffness_only"

# make random potential offsets hi
ϕx = 0
ϕy = 0
ϕz = 0

# QUASIPERIODIC
# p-wave 
Js = [2.74, 3.0, 3.24, 3.49, 3.74, 4.24, 5.24]
for J in Js
    V0 = 0
    V1 = -3
    ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="48:00:00")
end

# d-wave 
Js = [1.75, 2, 3.24, 4.49, 4.74, 4.99, 5.24, 5.99]
for J in Js
    V0 = 1
    V1 = -3
    ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="48:00:00")
end

# s-wave 
Js = collect(0:0.25:6)
push!(Js, 3.76, 4.76)
for J in Js
    V0 = -3
    V1 = 0
    ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="4:00:00")
end

# for J in Js
#     # QUASIPERIODIC
#     # p-wave 
#     V0 = 0
#     V1 = -3
#     ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
#     submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="48:00:00")

#     # d-wave 
#     V0 = 1
#     V1 = -3
#     ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
#     submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="48:00:00")

#     # s-wave 
#     V0 = -3
#     V1 = 0
#     ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=false)
#     submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="4:00:00")


#     # DISORDER 
#     # d-wave 
#     # V0 = 1
#     # V1 = -3
#     # ps = ModelParams(L=13, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
#     # submit_job(ps, filepath, @__DIR__, job_prefix, mem=800, time="6:00:00")
#     # ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
#     # submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="4:00:00")
#     # ps = ModelParams(L=7, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
#     # submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="10:00")

#     # # s-wave 
#     # V0 = -3
#     # V1 = 0
#     # ps = ModelParams(L=13, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
#     # submit_job(ps, filepath, @__DIR__, job_prefix, mem=800, time="6:00:00")
#     # ps = ModelParams(L=11, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
#     # submit_job(ps, filepath, @__DIR__, job_prefix, mem=512, time="4:00:00")
#     # ps = ModelParams(L=7, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=true)
#     # submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="10")
# end

