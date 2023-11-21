## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots, StatsBase
using Plots.PlotMeasures
include("../../src/model.jl")
include("utilities.jl")

# Parameters 
L = 11
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
periodic = true
disorder = true
savefigs = true
ndims = 3
nbins = 30

mkpath(joinpath(@__DIR__, "figures"))

if disorder
    pot = "disorder"
else
    pot = "QP"
end

# new data 
# dirname = "$(ndims)D_data_$pot"
# all_dfs = load_dfs(dirname=dirname)

# old data 
all_dfs = DataFrame(CSV.File("/Users/nicole/Dropbox/Grad/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/data_collection/IPR/3D_11L_data/IPR_data_11L.csv"))
df = all_dfs[(all_dfs.L.==L).&&(all_dfs.pot.==pot), :]
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

minE, maxE = -2 * ndims, 2 * ndims
binsize_E = 2 * maxE / nbins
binsize_ipr = 1 / nbins
gdf = groupby(df, [:J])
df_mean = DataFrame(J=[], pot=[], E=[], ipr_real_mean=[], ipr_k_mean=[], ipr_real_sem=[], ipr_k_sem=[])
for g in gdf
    iprs_r, iprs_k, Es = hcat(g.ipr_real...), hcat(g.ipr_k...), hcat(g.E...)
    # we need to rescale E by J 
    J = g.J[1]
    Es = Es ./ (1 + J / 2)
    iprs_r_avg, iprs_k_avg = mean(iprs_r, dims=2), mean(iprs_k, dims=2)
    iprs_r_sem, iprs_k_sem = sem_dims(iprs_r, dims=2), sem_dims(iprs_k, dims=2)

    Es_binned = zeros(size(Es)[2], nbins)
    for (Eᵢ, E) in enumerate(eachcol(Es))
        h = StatsBase.fit(Histogram, E, minE:binsize_E:maxE)
        Es_binned[Eᵢ, :] = h.weights ./ binsize_E
    end
    Es_binned_avg = mean(Es_binned, dims=1)
    Es_binned_avg = Es_binned_avg[1, :] # mean returns 2 dims still 

    #Es = StatsBase.fit(Histogram, E_avg[:, 1], minE:binsize_E:maxE).weights #./ binsize_E
    iprs_k_avg = bin_array_evenly(iprs_k_avg, nbins)
    iprs_r_avg = bin_array_evenly(iprs_r_avg, nbins)
    # Es = hist_counts(Es, nbins=nbins)
    # iprs_k_avg = hist_counts(iprs_k_avg, nbins=nbins)
    # iprs_r_avg = hist_counts(iprs_r_avg, nbins=nbins)

    dfi = DataFrame(J=[g.J[1]], pot=[g.pot[1]], E=[Es_binned_avg], ipr_real_mean=[iprs_r_avg], ipr_k_mean=[iprs_k_avg], ipr_real_sem=[iprs_r_sem], ipr_k_sem=[iprs_k_sem])
    append!(df_mean, dfi)
end

subdf = df_mean[(df_mean.pot.==pot), :]
dos, iprs_real, iprs_k, Js = hcat(subdf.E...), hcat(subdf.ipr_real_mean...), hcat(subdf.ipr_k_mean...), subdf.J

iprs = iprs_real
nx = length(Js)
ny = nbins
hval = copy(iprs)'
dos = dos'
E_edges = collect(minE:binsize_E:maxE)[1:end-1]

function get_colour(val; max_val)
    c = floor(Int, val / max_val * 100)
    if c == 0
        c = 1
    end
    # palett = palette([:purple, :red], 100)
    cmap = reverse(cgrad(:seismic, 200, categorical=true)[1:100])
    return cmap[c]
end

arr = zeros(4, 4)
p = heatmap(arr, alpha=0, c=:seismic, clims=(-1, 1))

plot!(p, grid=false, xlabel="J", ylabel="E/(1+J/2)", xticks=sort(Js)[1:3:end], margin=5mm, xaxis=:log10)
for (x, J) in Iterators.reverse(enumerate(Js)) # potential strength  
    for (y, E) in Iterators.reverse(enumerate(E_edges)) # eigenstates 
        dos_xy = dos[x, y]
        val_xy = hval[x, y]
        # phase 
        c = get_colour(val_xy, max_val=maximum(hval))
        # onsite term  
        scatter!(p, [J], [E], ms=abs(dos_xy) * 0.03, c=c, legend=:false)
    end
end

iprs = iprs_k
hval = copy(iprs)'

function get_colour(val; max_val)
    c = floor(Int, val / max_val * 100)
    if c == 0
        c = 1
    end
    # palett = palette([:purple, :blue], 100)
    cmap = cgrad(:seismic, 200, categorical=true)[100:end]
    return cmap[c]
end

for (x, J) in Iterators.reverse(enumerate(Js)) # potential strength  
    for (y, E) in Iterators.reverse(enumerate(E_edges)) # eigenstates 
        if J < 2
            dos_xy = dos[x, y]
            val_xy = hval[x, y]
            # phase 
            c = get_colour(val_xy, max_val=maximum(hval))
            # onsite term  
            scatter!(p, [J], [E], ms=abs(dos_xy) * 0.03, c=c, legend=:false)
        end
    end
end
xtcks = collect(0.2:0.9:20)
plot!(p, xticks=(xtcks, xtcks ./ 2))

if savefigs
    savefig(p, joinpath(@__DIR__, "figures", "$(L)L_$(ndims)D_IPR_$(pot).pdf"))
end
