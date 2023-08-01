
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames
using CurveFit

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")

function fit_χ(Ts, χs)
    # split up Ts into 5 bits 
    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χs = χs[sortidx]

    wl = ceil(Int, length(Ts) / 5) # window length 
    errs = []
    for idx in 1:(length(Ts)-wl)
        T = log10.(Ts[idx:idx+wl])
        χ = χs[idx:idx+wl]
        b, a = linear_fit(T, χ)
        χ̂ = a .* T .+ b
        err = norm(χ - χ̂)
        push!(errs, err)
    end

    minerr = argmin(errs)
    T = log10.(Ts[minerr:minerr+wl])
    χ = χs[minerr:minerr+wl]
    b, a = linear_fit(T, χ)
    χ̂ = a .* T .+ b

    return 10 .^ T, χ̂, a, b, errs[minerr]
end