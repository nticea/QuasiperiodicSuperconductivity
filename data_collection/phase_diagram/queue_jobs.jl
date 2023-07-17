## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include("../../src/model.jl")
include("../../src/submit_job.jl")
include("utilities.jl")

L = 17 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
ϕx = 0
ϕy = 0
periodic = 1

Js = collect(0:0.1:2)
V0s = collect(-1:0.1:1)
V1s = [-1.0, -0.5, 0]#collect(-2:0.5:0)
filepath = joinpath(@__DIR__, "collect_data.jl")
job_prefix = "phase_diagram"
dfs = load_dfs()

for J in Js
    for V0 in V0s
        for V1 in V1s
            if !already_computed(dfs, L=L, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, J=J)
                ps = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, V0=V0, V1=V1, J=J, periodic=periodic)
                submit_job(ps, filepath, @__DIR__, job_prefix)
            end
        end
    end
end