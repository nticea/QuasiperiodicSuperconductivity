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
include("../../src/BdG_dwave.jl")

## PARAMETERS 

L = 17
V0 = 1.5
V1 = -1

## SCRIPT 

df_BdG = load_dataframe("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_BdG_comparison/old_data/$(L)L$(V0)V0$(V1)V1_BdG_results.csv")
df_LGE = load_dataframe("/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/LGE_BdG_comparison/old_data/$(L)L$(V0)V0$(V1)V1_LGE_results.csv")

# function argmin_greater_than_val(arr; val)
#     posvals = findall(x -> x > val, arr)
#     idx = argmin(abs.(arr[posvals] .- val))
#     λmin = arr[posvals][idx]
#     findall(x -> x == λmin, arr)[1]
# end

# global hmaps = []
# Js = unique(df_LGE.J)

# for J in Js
#     print(J, "-")
#     # get the corresponding data
#     subdf_LGE = df_LGE[(df_LGE.J.==J).&(df_LGE.V0.==V0).&(df_LGE.V1.==V1), :]
#     subdf_BdG = df_BdG[(df_BdG.J.==J).&(df_BdG.V0.==V0).&(df_BdG.V1.==V1), :]

#     # get the LGE result 
#     λs = subdf_LGE.λ
#     if length(λs) > 0
#         idx = argmin_greater_than_val(λs, val=1)
#         Δ_LGE = subdf_LGE.Δ[idx]
#         evs = zeros(5, L, L)
#         for (n, i) in enumerate(1:(L*L):(5*L*L))
#             evi = Δ_LGE[i:(i+L*L-1)]
#             evs[n, :, :] = reshape(evi, L, L)
#         end
#         Δ_LGE = copy(evs)
#         λ = subdf_LGE.λ[idx]

#         # get the corresponding BdG result 
#         T = subdf_LGE.T[idx]
#         subdf_BdG = subdf_BdG[(subdf_BdG.T.==T), :]
#         @assert size(subdf_BdG)[1] == 1
#         Δ_BdG = subdf_BdG.Δ[1]
#         evs = zeros(5, L, L)
#         for (n, i) in enumerate(1:(L*L):(5*L*L))
#             evi = Δ_BdG[i:(i+L*L-1)]
#             evs[n, :, :] = reshape(evi, L, L)
#         end
#         Δ_BdG = copy(evs)

#         # plot both results 
#         θ = θ_to_π(subdf_LGE.θ[idx])
#         ϕx = θ_to_π(subdf_LGE.ϕx[idx])
#         ϕy = θ_to_π(subdf_LGE.ϕy[idx])

#         function colour_phase(x1::Int, x2::Int, x3::Int; all_evs, numpts::Int=10)
#             val = all_evs[x1, x2, x3]
#             if val < 0
#                 return "blue"
#             else
#                 return "red"
#             end
#         end

#         evs = Δ_LGE
#         p_LGE = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
#         if λ > 0
#             for x in 1:L
#                 for y in 1:L

#                     # bonds 
#                     plot!(p_LGE, [x, x - 1], [y, y], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
#                     plot!(p_LGE, [x, x], [y, y + 1], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
#                     plot!(p_LGE, [x, x + 1], [y, y], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
#                     plot!(p_LGE, [x, x], [y, y - 1], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)

#                     # onsite dot 
#                     scatter!(p_LGE, [x], [y], ms=100 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)
#                 end
#             end
#         end

#         xlabel!(p_LGE, "Site (x)")
#         ylabel!(p_LGE, "Site, (y)")
#         title!(p_LGE, "LGE sol'n T=$(round(T,digits=4)), λ=$(round(λ,digits=2)) \n Δ(J=$J, θ=$θ, ϕx=$ϕx, ϕy=$ϕy, \n V0=$V0, V1=$(round(V1,digits=2)))", fontsize=4)

#         evs = Δ_BdG .* (Δ_LGE[1, 1, 1] / Δ_BdG[1, 1, 1])
#         p_BdG = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
#         if λ > 0
#             for x in 1:L
#                 for y in 1:L

#                     # bonds 
#                     plot!(p_BdG, [x, x - 1], [y, y], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
#                     plot!(p_BdG, [x, x], [y, y + 1], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
#                     plot!(p_BdG, [x, x + 1], [y, y], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
#                     plot!(p_BdG, [x, x], [y, y - 1], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)

#                     # onsite dot 
#                     scatter!(p_BdG, [x], [y], ms=100 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)
#                 end
#             end
#         end

#         xlabel!(p_BdG, "Site (x)")
#         ylabel!(p_BdG, "Site, (y)")
#         title!(p_BdG, "BdG sol'n T=$(round(T,digits=4)), λ=$(round(λ,digits=2)) \n Δ(J=$J, θ=$θ, ϕx=$ϕx, ϕy=$ϕy, \n V0=$V0, V1=$(round(V1,digits=2)))", fontsize=4)

#         p = plot(p_LGE, p_BdG, layout=Plots.grid(1, 2, widths=[1 / 2, 1 / 2]), size=(1500, 1000), aspect_ratio=:equal)
#         push!(hmaps, p)

#         savefig(p, "comparison_$(J)J.png")
#     end
# end

Tc_BdG = find_Tc_BdG(df_BdG, interp_value=1e-6)
Tc_LGE = find_Tc_LGE(df_LGE, interp_value=1)

plot(Tc_BdG.J, Tc_BdG.Tc, label=nothing, c="blue")
scatter!(Tc_BdG.J, Tc_BdG.Tc, label="BdG soln", c="blue")
plot!(Tc_LGE.J, Tc_LGE.Tc, label=nothing, c="red")
scatter!(Tc_LGE.J, Tc_LGE.Tc, label="LGE soln", c="red")
title!("Comparing BdG and LGE solutions \n L=$L, V0=$V0, V1=$V1")

Js = unique(df_LGE.J)
p = plot()
for (i, J) in enumerate(Js[13])

    # Get the BdG data 
    dfsub_BdG = df_BdG[(df_BdG.J.==J).&(df_BdG.V0.==V0).&(df_BdG.V1.==V1), :]
    Δs, Ts = dfsub_BdG.Δ, dfsub_BdG.T

    max_Δ = []
    for Δ in Δs
        push!(max_Δ, maximum(Δ))
    end

    plot!(p, Ts, max_Δ, label=nothing, c="blue", yaxis=:log10)
    scatter!(p, Ts, max_Δ, label="J=$J", c="blue", yaxis=:log10)

    # Get the LGE data 
    dfsub_LGE = df_LGE[(df_LGE.J.==J).&(df_LGE.V0.==V0).&(df_LGE.V1.==V1), :]
    Tc = find_Tc_LGE(dfsub_LGE).Tc

    vline!(p, [Tc], label="LGE Tc", c="red")
end
title!(p, "BdG and LGE comparison")
xlabel!(p, "T")
ylabel!(p, "Maximum Δ")
plot(p)
