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
        if size(dfi)[1] == 1 && size(dfi)[2] == size(df_LGE)[2]
            dfi = convert_df_arrays(dfi, "Δ")
            dfi = convert_df_arrays(dfi, "K")
            dfi = convert_df_arrays(dfi, "Π")
            append!(df_LGE, dfi)
        end
    end
end

# extract only the parameters we are interested in 
df_LGE = df_LGE[(df_LGE.L.==L).&(df_LGE.θ.==θ).&(df_LGE.Q.==Q).&(df_LGE.ϕx.==ϕx).&(df_LGE.ϕy.==ϕy), :]
clims = (finite_minimum(df_LGE.T), finite_maximum(df_LGE.T))

## FIGURES ##
ps = []
grouped = groupby(df_LGE, [:V1])
V1s = []
for (i, df) in enumerate(grouped)
    mat = IndexMatrix(unique(df_LGE.V0), unique(df_LGE.J))

    V1 = df.V1[1]
    V0s, Js = sort(df.V0), sort(df.J)
    push!(V1s, V1)

    for V0 in V0s
        for J in Js
            dfi = df[(df.V0.==V0).&(df.J.==J).&(df.V1.==V1), :]
            if size(dfi)[1] == 1
                mat[V0, J] = dfi.T[1]
            else
                @error "Perhaps we need to average"
            end
        end
    end

    xtick_positions = collect(1:length(mat.cols))
    xtick_labels = [string(round(x, digits=2)) for x in mat.cols]
    ytick_positions = collect(1:length(mat.rows))
    ytick_labels = [string(round(y, digits=2)) for y in mat.rows]

    h = heatmap(mat.data, title="V1=$V1", xticks=(xtick_positions, xtick_labels), yticks=(ytick_positions, ytick_labels), clims=clims)
    xlabel!("J")
    ylabel!("V0")
    push!(ps, h)
end
sortidx = sortperm(V1s)
ps = ps[sortidx]

p = plot(ps..., layout=Plots.grid(2, 3,
        widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1000), plot_title="LGE Tc at θ=$(θ_to_π(θ))) for $L × $L lattice")
