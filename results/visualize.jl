## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using CSV
using DataFrames
using StatsPlots

L = 49
files = readdir(@__DIR__)

# nodenames = ["L", "J", "V0", "T", "λ"]
# df = DataFrame([name => [] for name in nodenames])
# for f in files
#     if startswith(f, "$(L)Nx$(L)Ny_results")
#         dfi = DataFrame(CSV.File(joinpath(@__DIR__, f)))
#         append!(df, dfi)
#     end
# end
loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/results/19Nx19Ny_results.csv"
df = DataFrame(CSV.File(loadpath))

Tc_df = find_Tc(df)

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