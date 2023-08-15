## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("../../src/submit_job.jl")

t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
ϕx = 0
ϕy = 0
ϕz = 0
periodic = 1
ndims = 3

Ls = [35, 30, 25, 20]
Js = collect(0:0.25:3)

filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "susceptibility"

for L in Ls
    for J in Js
        ps = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
        submit_job(ps, filepath, @__DIR__, job_prefix, mem=500)
    end
end