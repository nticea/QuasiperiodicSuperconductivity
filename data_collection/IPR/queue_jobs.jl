## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("submit_job.jl")

L = 23
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
ndims = 3

# L = 29, time = 6 hrs, mem = 350
# L = 27, time = 4 hrs, mem = 256
# L = 23, time = 1.5 hrs, mem = 100

Js = expspace(log10(1) - 1, log10(1) + 1, 20)

filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "IPR"

nrep = 40

for _ in 1:nrep
    # make random potential offsets 
    ϕx = 2π * rand()
    ϕy = 2π * rand()
    ϕz = 2π * rand()

    ## 3D 
    # QUASIPERIOIDIC!
    for J in Js
        ps = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=1, ndims=ndims, disorder=0)
        submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="30:00")
    end

    # DISORDER! 
    for J in Js
        ps = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=1, ndims=ndims, disorder=1)
        submit_job(ps, filepath, @__DIR__, job_prefix, mem=128, time="30:00")
    end
end
