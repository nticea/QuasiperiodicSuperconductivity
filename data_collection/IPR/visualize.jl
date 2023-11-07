## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")
include("../../src/utilities.jl")

# Parameters 
L = 11
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
periodic = true
ndims = 3

savepath = joinpath(@__DIR__, "data", "IPR_data_$(L)L.csv")
df = DataFrame(CSV.File(savepath))
df = convert_df_arrays(df, "ipr_real")
df = convert_df_arrays(df, "ipr_k")
df = convert_df_arrays(df, "E")

function sem_dims(arr; dims)
    @assert Base.ndims(arr) == 2
    res = []
    if dims == 1
        for a in eachcol(arr)
            push!(res, sem(a))
        end
    elseif dims == 2
        for a in eachrow(arr)
            push!(res, sem(a))
        end
    else
        println("idk man :(")
        return
    end
    return res
end

gdf = groupby(df, [:J, :pot])
df_mean = DataFrame(J=[], pot=[], ipr_real_mean=[], ipr_k_mean=[], ipr_real_sem=[], ipr_k_sem=[])
for g in gdf
    iprs_r, iprs_k = hcat(g.ipr_real...), hcat(g.ipr_k...)
    iprs_r_avg, iprs_k_avg = mean(iprs_r, dims=2), mean(iprs_k, dims=2)
    iprs_r_sem, iprs_k_sem = sem_dims(iprs_r, dims=2), sem_dims(iprs_k, dims=2)
    dfi = DataFrame(J=[g.J[1]], pot=[g.pot[1]], ipr_real_mean=[iprs_r_avg], ipr_k_mean=[iprs_k_avg], ipr_real_sem=[iprs_r_sem], ipr_k_sem=[iprs_k_sem])
    append!(df_mean, dfi)
end

