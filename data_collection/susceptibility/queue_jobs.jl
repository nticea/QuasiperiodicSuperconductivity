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
periodic = 1

Ls = [105]
Js = [0.1, 0.15, 0.17, 0.18, 0.19, 0.21, 0.22, 0.23, 0.25, 0.3]#collect(0:0.2:3)

filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "susceptibility"

for L in Ls
    for J in Js
        ps = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, J=J, periodic=periodic)
        submit_job(ps, filepath, @__DIR__, job_prefix, mem=256)
    end
end