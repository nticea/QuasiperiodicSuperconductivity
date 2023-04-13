## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
include("../src/BdG.jl")
include("../src/results.jl")

using Plots
using CSV
using DataFrames
using StatsPlots

loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/BdG_results/9Nx9Ny_results.csv"
df = DataFrame(CSV.File(loadpath))

## PARAMETERS ##
L = 9 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
pairing_symmetry = "s-wave"
tol = 1e-6

# read files
# files = readdir(@__DIR__)
# nodenames = ["L", "J", "V0", "T", "λ"]
# df = DataFrame([name => [] for name in nodenames])
# for f in files
#     if startswith(f, "$(L)Nx$(L)Ny_results")
#         dfi = DataFrame(CSV.File(joinpath(@__DIR__, f)))
#         append!(df, dfi)
#     end
# end

@assert length(unique(df.L)) == 1
@assert df.L[1] == L

ΔE = finite_size_gap(L=L, t=t, J=0, Q=Q, μ=μ)
fsgap = maximum(ΔE)

Tc_df = find_Tc(df, interp_value=fsgap)

# Now, do plotting 
@df Tc_df plot(
    :V0,
    :Tc,
    group=:J,
    m=(0.75, [:+ :h :star7 :circle], 5),
    xaxis=:log10, yaxis=:log10,
)
title!("Transition temperature")
xlabel!("V")
ylabel!("Tc")