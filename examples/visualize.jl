## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using CSV
using DataFrames
using StatsPlots

# loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/29Nx29Ny_results.csv"
# loadpath2 = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/29Nx29Ny_results_copy.csv"
# df = DataFrame(CSV.File(loadpath))
# df_old = DataFrame(CSV.File(loadpath2))

# append!(df, df_old)

loadpath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/results/49Nx49Ny_results.csv"
df = DataFrame(CSV.File(loadpath))

# Get the unique Js, V0s from the dataframe 
Js, V0s = unique(df.J), unique(df.V0)
L = df.L[1]

# Make the new dataframe 
nodenames = ["L", "J", "V0", "Tc"]
Tc_df = DataFrame([name => [] for name in nodenames])

# For each unique (J,V0) pair... extract the corresponding data across all Ts 
for J in Js
    for V0 in V0s

        # Extract the corresponding data across all Ts
        dfsub = df[(df.J.==J).&(df.V0.==V0), :]
        λs, Ts = dfsub.λ, dfsub.T

        # Compute the interpolated Tc for this (J,V0) pair
        knots = reverse(λs)
        if length(knots) > 1
            Interpolations.deduplicate_knots!(knots, move_knots=true)
            try
                interp_linear = linear_interpolation(knots, reverse(Ts))
                Tc = interp_linear(1)
                # Put it into a new dataframe indexed by (J,V0,Tc)
                df2 = DataFrame(L=[L], Tc=[Tc], J=[J], V0=[V0])
                append!(Tc_df, df2)
            catch e
                Tc = NaN
                # Put it into a new dataframe indexed by (J,V0,Tc)
                df2 = DataFrame(L=[L], Tc=[Tc], J=[J], V0=[V0])
                append!(Tc_df, df2)
            end
        end
    end
end

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