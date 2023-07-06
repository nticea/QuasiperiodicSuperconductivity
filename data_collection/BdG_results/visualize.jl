## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")

## PARAMETERS ##
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
periodic = true
T = 0
tol = 1e-15

# load in the data 
savepath_BdG = joinpath(@__DIR__, "BdG_results.csv")
df_BdG = load_dataframe(savepath_BdG)
