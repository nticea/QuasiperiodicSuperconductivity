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
savepath_LGE = joinpath(@__DIR__, "LGE_Tc.csv")
df = load_dataframe(savepath_BdG)
Tc_df = load_dataframe(savepath_LGE)

# Get the unique Js, V1s from the dataframe 
Ls, θs, Js, V0s, V1s = unique(df.L), unique(df.θ), unique(df.J), unique(df.V0), unique(df.V1)

function convert_df_arrays(df::DataFrame, col_name::String, delims=r"[,; ]")
    all_arrs = []
    for (i, arr) in enumerate(df[!, Symbol(col_name)])
        # Split the string by multiple delimiters
        result = split(chop(arr; head=1, tail=1), delims)
        new_array = result[result.!=""]
        arrs_new = parse.(Float64, new_array)
        push!(all_arrs, arrs_new)
    end
    dfcut = copy(df)#df[:, collect(1:size(df)[2]-1)]
    dfcut[!, Symbol(col_name)] = all_arrs

    return dfcut
end

df = convert_df_arrays(df, "K")
df = convert_df_arrays(dfnew, "Π")

# Make the new dataframe 
nodenames = ["L", "θ", "J", "V0", "V1", "K", "Π", "Tc"]
df_avg = DataFrame([name => [] for name in nodenames])

# For each unique (J,V0) pair... extract the corresponding data across all Ts 
for L in Ls
    for θ in θs
        for J in Js
            for V0 in V0s
                for V1 in V1s
                    # Extract the corresponding data
                    dfsub = df[(df.L.==L).&(df.θ.==θ).&(df.J.==J).&(df.V0.==V0).&(df.V1.==V1), :]
                    dfsub_Tc = Tc_df[(Tc_df.L.==L).&(Tc_df.θ.==θ).&(Tc_df.J.==J).&(Tc_df.V0.==V0).&(Tc_df.V1.==V1), :]
                    Ks, Πs = dfsub.K, dfsub.Π
                    Ks, Πs = hcat(Ks...), hcat(Πs...)
                    Ts = dfsub_Tc.T

                    # take the mean 
                    Ks, Πs = mean(Ks, dims=2)[:, 1], mean(Πs, dims=2)[:, 1]
                    Ts = mean(Ts)

                    df2 = DataFrame(L=[L], θ=[θ], J=[J], V0=[V0], V1=[V1], K=[Ks], Π=[Πs], Tc=[Ts])
                    append!(df_avg, df2)
                end
            end
        end
    end
end

# now it is time to plot 
Js, Ks, Πs, Tcs = df_avg.J, df_avg.K, df_avg.Π, df_avg.Tc
Ds = (-Ks + Πs)
Ds = hcat(Ds...)

L, θ, V0, V1 = df_avg.L[1], df_avg.θ[1], df_avg.V0[1], df_avg.V1[1]

p = plot()
dirs = ["x̂", "ŷ", "-x̂", "-ŷ"]
cmap = cgrad(:matter, 4, categorical=true)

for i in 1:4
    plot!(p, Js, Ds[i, :], label=nothing, c=cmap[i])
    scatter!(p, Js, Ds[i, :], label=dirs[i], c=cmap[i])
end
plot!(Js, Tcs, color="blue", label=nothing)
scatter!(Js, Tcs, color="blue", label="LGE soln")
title!(p, "Superfluid stiffness for \n Δ(V0=$V0, V1=$V1, θ=$(θ_to_π(θ))) \n $L × $L lattice at T=0")
xlabel!("J")
ylabel!("Tc")