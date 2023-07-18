## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using CSV
using DataFrames
using StatsPlots

include("../../src/model.jl")
include("../../src/BdG_dwave.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")

## PARAMETERS ## 

L = 17 # the full system is L × L 
Q = (√5 - 1) / 2
θ = π / 7
V0_swave = -2.3
V1_swave = 0
V0_dwave = 1
V1_dwave = -1.5
ϕx = 0
ϕy = 0
savefigs = true

# read files 
files = readdir(joinpath(@__DIR__, "data"))
df_BdG = DataFrame(L=Int64[], J=Float64[], Q=Float64[], θ=Float64[],
    ϕx=Float64[], ϕy=Float64[], V0=Float64[], V1=Float64[],
    T=Float64[], λ=Float64[], Δ=[], K=[], Π=[])
df_LGE = copy(df_BdG)
for f in files
    if endswith(f, ".csv")
        dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
        # process all of the arrays 
        dfi = convert_df_arrays(dfi, "Δ")
        dfi = convert_df_arrays(dfi, "K")
        dfi = convert_df_arrays(dfi, "Π")
        if contains(f, "BdG")
            append!(df_BdG, dfi)
        elseif contains(f, "LGE")
            append!(df_LGE, dfi)
        end
    end
end

# extract only the parameters we are interested in 
df_LGE_swave = df_LGE[(df_LGE.L.==L).&(df_LGE.θ.==θ).&(df_LGE.Q.==Q).&(df_LGE.V0.==V0_swave).&(df_LGE.V1.==V1_swave), :]
df_LGE_dwave = df_LGE[(df_LGE.L.==L).&(df_LGE.θ.==θ).&(df_LGE.Q.==Q).&(df_LGE.V0.==V0_dwave).&(df_LGE.V1.==V1_dwave), :]
df_BdG_swave = df_BdG[(df_BdG.L.==L).&(df_BdG.θ.==θ).&(df_BdG.Q.==Q).&(df_BdG.V0.==V0_swave).&(df_BdG.V1.==V1_swave), :]
df_BdG_dwave = df_BdG[(df_BdG.L.==L).&(df_BdG.θ.==θ).&(df_BdG.Q.==Q).&(df_BdG.V0.==V0_dwave).&(df_BdG.V1.==V1_dwave), :]

## FIGURES ##
# Comparing Tc from LGE and from BdG 
df_LGE_Tc_swave = df_LGE_swave[df_LGE_swave.T.>0, :]
df_LGE_Tc_dwave = df_LGE_dwave[df_LGE_dwave.T.>0, :]
df_BdG_Tc_swave = df_BdG_swave[df_BdG_swave.T.>0, :]
df_BdG_Tc_dwave = df_BdG_dwave[df_BdG_dwave.T.>0, :]

p1 = plot(size=(800, 600))
Js, Tcs = df_LGE_Tc_swave.J, df_LGE_Tc_swave.T
sortidx = sortperm(Js)
Js, Tcs = Js[sortidx], Tcs[sortidx]
plot!(p1, Js, Tcs, color="red", label=nothing)
scatter!(p1, Js, Tcs, color="red", label="s-wave Tc (V0=$V0_swave)", ylabel="Tc")

@assert 1 == 0

Js = sort(unique(df_LGE_Tc_dwave.J))
Tcs = zeros(length(Js))
for (j, J) in enumerate(Js)
    dfsub = df_LGE_Tc_dwave[(df_LGE_Tc_dwave.J.==J), :]
    Ts = dfsub.T
    Tcs[j] = mean(Ts)
end
plot!(p1, Js, Tcs, color="blue", label=nothing)
scatter!(p1, Js, Tcs, color="blue", label="d-wave Tc (V0=$V0_dwave, V1=$V1_dwave)", ylabel="Tc",)

# make some adjustments
plot!(p1, legend=:bottomright, ylims=(0, 0.37))
xlabel!(p1, "J")
title!(p1, "Comparing d-wave and s-wave LGE Tc on $L × $L lattice")

# Superfluid stiffness vs Tc 
# df_LGE_0 = df_LGE[df_LGE.T.==0, :]
# df_BdG_0 = df_BdG[df_BdG.T.==0, :]

# Js = sort(unique(df_BdG_0.J))
# Ds = zeros(length(Js), 4)
# for (j, J) in enumerate(Js)
#     dfsub = df_BdG_0[(df_BdG_0.L.==L).&(df_BdG_0.J.==J).&(df_BdG_0.θ.==θ).&(df_BdG_0.Q.==Q).&(df_BdG_0.V0.==V0).&(df_BdG_0.V1.==V1), :]
#     Ks, Πs = hcat(dfsub.K...), hcat(dfsub.Π...)
#     Ks, Πs = mean(Ks, dims=2), mean(Πs, dims=2)
#     Ds[j, :] = -Ks + Πs
# end

# dirs = ["Dₛ/π (x̂)", "Dₛ/π (ŷ)"]
# cmap = ["green", "orange"]
# for i in 1:2
#     plot!(p1, Js, Ds[:, i], label=nothing, c=cmap[i], secondary=true)
#     scatter!(p1, Js, Ds[:, i], label=dirs[i], c=cmap[i], secondary=true)
# end

if savefigs
    savefig(p1, joinpath(@__DIR__, "figures", "$(L)L_$(V0_swave)V0swave_$(V0_dwave)V0dwave_$(V1_dwave)V1dwave_Tc1.pdf"))
end

