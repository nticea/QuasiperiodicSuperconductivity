## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots, StatsPlots
using CSV, DataFrames
using CurveFit

include("../../src/model.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")
include("utilities.jl")

## PARAMETERS ## 

Ls = [89] # the full system is L × L 
Q = (√5 - 1) / 2
θ = π / 7
savefigs = true
figpath = mkpath(joinpath(@__DIR__, "figures"))

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
files = readdir(joinpath(@__DIR__, "data"))
df = DataFrame(L=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[],
    T=[], χ=[])
for f in files
    if endswith(f, ".csv")
        dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
        dfi = convert_df_arrays(dfi, "χ")
        append!(df, dfi)
    end
end

if length(Ls) == 1
    cmap = ["red"]
else
    cmap = cgrad(:matter, length(Ls), categorical=true)
end

px, pos = plot(), plot()
for (l, L) in enumerate(Ls)
    # extract only the parameters we are interested in 
    dfL = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q), :]
    Js = sort(unique(dfL.J))

    for (j, J) in enumerate(Js)
        dfJ = dfL[(dfL.J.==J), :]
        Ts = dfJ.T
        χs = dfJ.χ

        #sortidx = sortperm(Ts)
        #Ts = Ts[sortidx]
        #χx = χs[sortidx]

        if length(Ts) > 0

            # on-site
            χswave = [a[9] for a in χs]

            # make the d-wave components 
            xx = [a[1] for a in χs]
            yy = [a[5] for a in χs]
            xy = [a[2] for a in χs]
            yx = [a[4] for a in χs]
            χdwave = 1 / 4 * (xx + yy - xy - yx)

            #plot!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10)
            scatter!(px, Ts, χdwave, color=cmap[l], label="J=$J, L=$L", xaxis=:log10)
            # now plot the error
            T̂, χ̂, a, b, err = fit_χ(Ts, χdwave)
            #plot!(px, T̂, χ̂, color="red", xaxis=:log10, label="Error: $(round(err,digits=2))")

            #plot!(pos, Ts, χswave, color=cmap[j], label=nothing, xaxis=:log10)
            scatter!(pos, Ts, χswave, color=cmap[l], label="J=$J, L=$L", xaxis=:log10)
            # now plot the error
            T̂, χ̂, a, b, err = fit_χ(Ts, χswave)
            #plot!(pos, T̂, χ̂, color="red", xaxis=:log10, label="Error: $(round(err,digits=2))")

            title!(px, "Susceptibility (d-wave)")
            xlabel!(px, "T")
            ylabel!(px, "χ")
            title!(pos, "Susceptibility (s-wave)")
            xlabel!(pos, "T")
            ylabel!(pos, "χ")
        end
    end
end

p = plot(px, pos, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility")
if savefigs
    savefig(p, joinpath(figpath, "susceptibility.pdf"))
end