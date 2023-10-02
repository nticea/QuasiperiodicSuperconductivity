## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("submit_job.jl")

t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0.5
θ = π / 7
V0 = 0
V1 = 0
ϕx = 0
ϕy = 0
ϕz = 0
periodic = 1 # periodic 
ndims = 3
disorder = 0 # add disorder, just for this run! 

# L = 29, time = 6 hrs, mem = 350
# L = 27, time = 4 hrs, mem = 256
# L = 23, time = 1.5 hrs, mem = 100

Js = collect(0:0.1:3)
Ts = expspace(-3, 1, 30) # temperature 

filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "susceptibility"

nrep = 10

for _ in 1:nrep
    for J in Js
        for T in Ts
            ps = ModelParams(L=17, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=0)
            submit_job(ps, filepath, @__DIR__, job_prefix, mem=80, kwargs="$T", time="45:00")
        end
    end
end

for _ in 1:nrep
    for J in Js
        for T in Ts
            ps = ModelParams(L=17, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=1)
            submit_job(ps, filepath, @__DIR__, job_prefix, mem=80, kwargs="$T", time="45:00")
        end
    end
end


# for J in Js
#     for T in Ts
#         ps = ModelParams(L=15, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
#         submit_job(ps, filepath, @__DIR__, job_prefix, mem=75, kwargs="$T", time="30:00")
#     end
# end

# for J in Js
#     for T in Ts
#         ps = ModelParams(L=20, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
#         submit_job(ps, filepath, @__DIR__, job_prefix, mem=100, kwargs="$T", time="1:00:00")
#     end
# end

# for J in Js
#     for T in Ts
#         ps = ModelParams(L=25, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
#         submit_job(ps, filepath, @__DIR__, job_prefix, mem=255, kwargs="$T", time="5:00:00")
#     end
# end

# for J in Js
#     for T in Ts
#         ps = ModelParams(L=30, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
#         submit_job(ps, filepath, @__DIR__, job_prefix, mem=350, kwargs="$T", time="12:00:00")
#     end
# end
