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
using StatsBase

# loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/BdG_results/20Nx20Ny_results.csv"
# df = DataFrame(CSV.File(loadpath))

## PARAMETERS ##
L = 55 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
pairing_symmetry = "s-wave"
tol = 1e-4
periodic = true

# read files
files = readdir(@__DIR__)
nodenames = ["L", "J", "V0", "T", "λ"]
df = DataFrame([name => [] for name in nodenames])
for f in files
    if startswith(f, "$(L)Nx$(L)Ny_results")
        dfi = DataFrame(CSV.File(joinpath(@__DIR__, f)))
        append!(df, dfi)
    end
end

@assert length(unique(df.L)) == 1
@assert df.L[1] == L

ΔE = finite_size_gap(L=L, t=t, Q=Q, μ=μ, periodic=periodic)
fsgap = maximum(ΔE)

J = 3
dfJ = df[(df.J.==J), :]
V0s = unique(dfJ.V0)
cs = palette([:purple, :black], length(V0s))

global p = plot()
for (i, V0) in enumerate(V0s)
    dfsub = df[(df.J.==J).&(df.V0.==V0), :]
    Ts = dfsub.T
    Δs = dfsub.λ
    global p = plot!(p, Ts, Δs, xscale=:log10, yscale=:log10, label=nothing, c=cs[i])
    global p = scatter!(p, Ts, Δs, xscale=:log10, yscale=:log10, label="V=$(round(V0,digits=3))", legend=:bottomleft, c=cs[i])
end
plot(p)
hline!([fsgap], ls=:dash, c="red", label="Finite size gap")
title!("|Δ| for J=$J")
xlabel!("Temperature")

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