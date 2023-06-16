## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/results.jl")
include("../../src/meanfield.jl")

## PARAMETERS 

L = 17
# symmetry = "dwave"

# ## SCRIPT 

# loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_sweep_J/$(L)Nx$(L)Ny_$(symmetry)_results.csv"
# df = load_dataframe(loadpath)
# df = sort(df, :J)

# V0 = df.V0[1]
# V1 = df.V1[1]

# # for every unique V1, find the Δ with the λ closest to 0 
# Js = unique(df.J)
# global hmaps = []
# for J in Js
#     # get the corresponding data
#     subdf = df[(df.J.==J), :]
#     λs = subdf.λ
#     idx = argmin(abs.(λs .- 1))
#     Ts = subdf.T

#     @show J, V0, V1
#     @show Ts
#     @show λs

#     println("")
#     push!(hmaps, plot_LGE_Δ(subdf; idx=idx))
# end
# p = plot(hmaps[1:20]..., layout=Plots.grid(4, 5, widths=[1 / 5 for _ in 1:5]), size=(2000, 1500), aspect_ratio=:equal)

symmetry = "swave"
loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_sweep_J/$(L)Nx$(L)Ny_$(symmetry)_results.csv"

files = readdir(@__DIR__)
nodenames = ["L", "J", "V0", "V1", "Tc"]
global Tc_dfs = DataFrame([name => [] for name in nodenames])
for f in files
    if startswith(f, "$(L)Nx$(L)Ny_$(symmetry)") && endswith(f, "csv")
        @show f
        df = load_dataframe(joinpath(@__DIR__, f))
        Tc_df = find_Tc(df, interp_value=1)
        Tc_df = sort(Tc_df, :J)
        append!(Tc_dfs, Tc_df)
    end
end

plot(Tc_dfs.J, Tc_dfs.Tc, label=nothing, color="red")
scatter!(Tc_dfs.J, Tc_dfs.Tc, label=symmetry, color="red")
title!("Transition temperature")
xlabel!("J")
ylabel!("Tc")

symmetry = "dwave"
loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_sweep_J/$(L)Nx$(L)Ny_$(symmetry)_results.csv"

files = readdir(@__DIR__)
nodenames = ["L", "J", "V0", "V1", "Tc"]
global Tc_dfs = DataFrame([name => [] for name in nodenames])
for f in files
    if startswith(f, "$(L)Nx$(L)Ny_$(symmetry)") && endswith(f, "csv")
        @show f
        df = load_dataframe(joinpath(@__DIR__, f))
        Tc_df = find_Tc(df, interp_value=1)
        Tc_df = sort(Tc_df, :J)
        append!(Tc_dfs, Tc_df)
    end
end

plot!(Tc_dfs.J, Tc_dfs.Tc, label=nothing, color="blue")
scatter!(Tc_dfs.J, Tc_dfs.Tc, label=symmetry, color="blue")
title!("Transition temperature")
xlabel!("J")
ylabel!("Tc")