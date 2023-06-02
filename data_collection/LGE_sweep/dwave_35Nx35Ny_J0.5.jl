## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/model.jl")
include("../../src/meanfield.jl")

## PARAMETERS ##
L = 17 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
J = 0.5
V1 = -1
periodic = true
pairing_symmetry = "d-wave"

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_$(J)J_dwave_results.csv")

# load in the dataframe, if it exists. If not, make a new one
df = load_dataframe(savepath)

function LGE_sweep(Ts; α::Real)
    ## RUNNING THE CODE ## 
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))

    # Compute V1, holding the ratio α≡V0/V1 fixed 
    V0 = V1 * α

    for i in iter # iterate through all temperatures
        T = Ts[i]
        # check whether this particular (J,T,V0) combo has been already computed 
        if !already_calculated(df; L=L, J=J, θ=θ, V0=V0, V1=V1, T=T)

            # calculate M at a given (J,T,V0, V1)
            @time λ, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, μ=μ, V0=V0, V1=V1, periodic=periodic, symmetry=pairing_symmetry)

            update_results!(df; L=L, λ=λ, J=J, θ=θ, V0=V0, V1=V1, T=T, Δ=Δ)
            CSV.write(savepath, df)
            flush(stdout)
        end
    end
end

αs = range(0.5, stop=1.3, length=9)
Ts = [0.09, 0.13, 0.16, 0.19]

for α in αs
    LGE_sweep(Ts, α=α)
end

loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_sweep/$(L)Nx$(L)Ny_$(J)J_dwave_results.csv"
df = load_dataframe(loadpath)

# for every unique V1, find the Δ with the λ closest to 0 
V0s = unique(df.V0)
global hmaps = []
for V0 in V0s
    # get the corresponding data
    subdf = df[(df.V0.==V0), :]
    λs = subdf.λ
    idx = argmin(abs.(abs.(λs) .- 1))
    @show λs
    @show idx
    println("")
    push!(hmaps, plot_LGE_Δ(subdf; idx=idx))
end
p = plot(hmaps..., layout=Plots.grid(3, 3, widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1500), aspect_ratio=:equal)