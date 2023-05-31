## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames
using StatsPlots

include("../../src/model.jl")
include("../../src/meanfield.jl")

loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_sweep/17Nx17Ny_0J_dwave_results.csv"

df = load_dataframe(loadpath)

# for every unique V1, find the Δ with the λ closest to 0 
V0s = unique(df.V0)
global hmaps = []
for V0 in V0s
    # get the corresponding data
    subdf = df[(df.V0.==V0), :]
    λs = subdf.λ
    idx = argmin(abs.(λs .- 1))
    @show λs
    @show idx
    println("")
    push!(hmaps, plot_LGE_Δ(subdf; idx=idx))
end
p = plot(hmaps..., layout=Plots.grid(3, 3, widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1500), aspect_ratio=:equal)

@assert 1 == 0
# Now, plot the Tcs 

files = readdir(@__DIR__)
nodenames = ["L", "J", "V0", "V1", "Tc"]
global Tc_dfs = DataFrame([name => [] for name in nodenames])
for f in files
    if startswith(f, "35Nx35Ny_") && endswith(f, "csv")
        @show f
        df = load_dataframe(joinpath(@__DIR__, f))
        Tc_df = find_Tc(df, interp_value=1)
        Tc_df = sort(Tc_df, :V1)
        append!(Tc_dfs, Tc_df)
    end
end

# Now, do plotting 
@df Tc_dfs plot(
    :V1,
    :Tc,
    group=:J,
    m=(0.75, [:+ :h :star7 :circle], 5),
    xaxis=:log10, yaxis=:log10,
)
title!("Transition temperature")
xlabel!("V1")
ylabel!("Tc")