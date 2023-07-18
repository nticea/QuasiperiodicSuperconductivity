## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

L = 17 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = -2.3
V1 = 0

Js = collect(0:0.1:4)
ϕxs, ϕys = [0], [0] #LinRange(0, 3, 3), LinRange(0, 3, 3)
filepath = joinpath(@__DIR__, "collect_data_fast.jl")

for ϕx in ϕxs
    for ϕy in ϕys
        for J in Js
            # silly 
            periodic = 1

            # make the model parameters
            empty!(ARGS)
            push!(ARGS, "$L")
            push!(ARGS, "$t")
            push!(ARGS, "$Q")
            push!(ARGS, "$μ")
            push!(ARGS, "$θ")
            push!(ARGS, "$ϕx")
            push!(ARGS, "$ϕy")
            push!(ARGS, "$V0")
            push!(ARGS, "$V1")
            push!(ARGS, "$J")
            push!(ARGS, "$periodic")

            include(filepath)
        end
    end
end