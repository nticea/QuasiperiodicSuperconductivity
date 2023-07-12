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
ϕx = 0
ϕy = 0

# read files 
files = readdir(joinpath(@__DIR__, "data"))
df_LGE = DataFrame(L=Int64[], J=Float64[], Q=Float64[], θ=Float64[],
    ϕx=Float64[], ϕy=Float64[], V0=Float64[], V1=Float64[],
    T=Float64[], λ=Float64[], Δ=[], K=[], Π=[])
for f in files
    if endswith(f, ".csv") && contains(f, "LGE")
        dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
        # process all of the arrays 
        dfi = convert_df_arrays(dfi, "Δ")
        dfi = convert_df_arrays(dfi, "K")
        dfi = convert_df_arrays(dfi, "Π")
        append!(df_LGE, dfi)
    end
end

# extract only the parameters we are interested in 
df_LGE = df_LGE[(df_LGE.L.==L).&(df_LGE.θ.==θ).&(df_LGE.Q.==Q).&(df_LGE.ϕx.==ϕx).&(df_LGE.ϕy.==ϕy), :]

## FIGURES ##
ps = []
grouped = groupby(df_LGE, [:V1])
for (i, df) in enumerate(grouped)
    mat = IndexMatrix(unique(df_LGE.J), unique(df_LGE.V0))

    V1 = df.V1[1]
    V0s, Js = sort(df.V0), sort(df.J)

    for V0 in V0s
        for J in Js
            dfi = df[(df.V0.==V0).&(df.J.==J).&(df.V1.==V1), :]
            if size(dfi)[1] == 1
                mat[J, V0] = dfi.T[1]
            else
                @error "Perhaps we need to average"
            end
        end
    end
    h = heatmap(mat.data)
    push!(ps, h)
end

# for (i, df) in enumerate(grouped)
#     ϕx, ϕy = unique(df.ϕx)[1], unique(df.ϕy)[1]
#     Js, Tcs = df.J, df.T
#     sortidx = sortperm(Js)
#     Js, Tcs = Js[sortidx], Tcs[sortidx]
#     plot!(p1, Js, Tcs, color=cmap[i], label=nothing)
#     scatter!(p1, Js, Tcs, color=cmap[i], label="(ϕx,ϕy)=($(round(ϕx, digits=2)),$(round(ϕy, digits=2)))", ylabel="Tc",)
# end
# xlabel!(p1, "J")
# title!(p1, "Tc for V0=$V0, V1=$V1, θ=$(θ_to_π(θ))\n on $L × $L lattice")
