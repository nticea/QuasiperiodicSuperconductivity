## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using ProgressBars

## PARAMETERS ##
L = 49 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
pairing_symmetry = "s-wave"

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_results.npy")

# J, V0, T
J = 2
V0 = 1
T = 1e-2
λ = @time λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
@show λ
flush(stdout)
println("Repeating due to compilation time")
λ = @time λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
flush(stdout)

results = Results(L, λ, J, V0, T)
save_structs(results, savepath)