## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("../../src/submit_job.jl")

L = 23 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
periodic = 1

Js = collect(0:0.1:3)

filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "susceptibility"

for J in Js
    ps = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, J=J, periodic=periodic)
    submit_job(ps, filepath, @__DIR__, job_prefix, mem=10)
end