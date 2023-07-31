
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "../.."))
    include("../../src/meanfield.jl")

    t = 1 # hopping 
    Q = (√5 - 1) / 2
    μ = 1e-8
    θ = π / 7
    V0 = 0
    V1 = 0
    ϕx = 0
    ϕy = 0
    periodic = true
    J = 0
    L = 35
end

@everywhere function _uniform_susceptibility(T; L, t, J,
    Q, θ, ϕx, ϕy, μ, periodic)

    uniform_susceptibility(T, L=L, t=t, J=J,
        Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic)
end

Ts = expspace(-4, 1, 5) # temperature 

# export the function
@sync @distributed for T in Ts
    χ = @time _uniform_susceptibility(T, L=L, t=t, J=J,
        Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, periodic=periodic)
end


