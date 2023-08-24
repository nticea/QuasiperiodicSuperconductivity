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

L = 23 # the full system is L × L 
ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
savefigs = true
figpath = mkpath(joinpath(@__DIR__, "figures"))
T_cutoff = 1e-2

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
files = readdir(joinpath(@__DIR__, "data"))
df = load_dfs()

px, pos = plot(margin=10Plots.mm), plot(margin=10Plots.mm)
# extract only the parameters we are interested in 
dfL = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q).&(df.ndims.==ndims), :]
Js = sort(unique(dfL.J))
cmap = cgrad(:matter, length(Js), categorical=true)

for (j, J) in enumerate(Js)
    dfJ = dfL[(dfL.J.==J), :]
    Ts = dfJ.T
    χs = dfJ.χ

    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χs = χs[sortidx]
    χs = [reshape(χ, 4, 4) for χ in χs]

    # consider only the linear regime 
    # χs = χs[Ts.>=T_cutoff]
    # Ts = Ts[Ts.>=T_cutoff]

    if length(Ts) > 0

        # on-site
        χswave = [χ[1, 1] for χ in χs]

        # make the d-wave components 
        xx, yy = [χ[2, 2] for χ in χs], [χ[3, 3] for χ in χs]
        xy, yx = [χ[2, 3] for χ in χs], [χ[3, 2] for χ in χs]
        if ndims == 2
            χdwave = xx + yy - xy - yx
        elseif ndims == 3
            zz = [χ[4, 4] for χ in χs]
            xz, zx = [χ[2, 4] for χ in χs], [χ[4, 2] for χ in χs]
            yz, zy = [χ[3, 4] for χ in χs], [χ[4, 3] for χ in χs]
            χdwave = xx + yy + zz - xy - yx - xz - zx - yz - zy
        else
            println("sorry")
            χdwave = nothing
        end

        plot!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10, yaxis=:log10)
        scatter!(px, Ts, χdwave, color=cmap[j], label="J=$J, L=$L", xaxis=:log10, yaxis=:log10)

        plot!(pos, Ts, χswave, color=cmap[j], label=nothing, xaxis=:log10, yaxis=:log10)
        scatter!(pos, Ts, χswave, color=cmap[j], label="J=$J, L=$L", xaxis=:log10, yaxis=:log10)

        title!(px, "d-wave")
        xlabel!(px, "T")
        ylabel!(px, "χ")
        title!(pos, "s-wave")
        xlabel!(pos, "T")
        ylabel!(pos, "χ")
    end
end

p1 = plot(px, pos, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility for $size_str lattice with Q=$(round(Q, digits=3)), θ=$(θ_to_π(θ))")
if savefigs
    savefig(p, joinpath(figpath, "susceptibility_J_$(L)L.pdf"))
end

# also, plot the value of each J at various temperatures 
Ts = sort(unique(dfL.T))
Ts = Ts[Ts.>=T_cutoff]
cmap = cgrad(:viridis, length(Ts), categorical=true)

ptemp_s = plot(margin=10Plots.mm)
ptemp_d = plot(margin=10Plots.mm)

for (i, T) in enumerate(Ts)
    dfT = dfL[(dfL.T.==T), :]
    Js = dfT.J
    χs = dfT.χ

    sortidx = sortperm(Js)
    Js = Js[sortidx]
    χs = χs[sortidx]
    χs = [reshape(χ, 4, 4) for χ in χs]

    # on-site
    χswave = [χ[1, 1] for χ in χs]

    # make the d-wave components 
    xx, yy = [χ[2, 2] for χ in χs], [χ[3, 3] for χ in χs]
    xy, yx = [χ[2, 3] for χ in χs], [χ[3, 2] for χ in χs]
    if ndims == 2
        χdwave = xx + yy - xy - yx
    elseif ndims == 3
        zz = [χ[4, 4] for χ in χs]
        xz, zx = [χ[2, 4] for χ in χs], [χ[4, 2] for χ in χs]
        yz, zy = [χ[3, 4] for χ in χs], [χ[4, 3] for χ in χs]
        χdwave = xx + yy + zz - xy - yx - xz - zx - yz - zy
    else
        println("sorry")
        χdwave = nothing
    end

    plot!(ptemp_s, Js, χswave, label=nothing, c=cmap[i])
    scatter!(ptemp_s, Js, χswave, c=cmap[i], label=nothing)

    plot!(ptemp_d, Js, χdwave, label=nothing, c=cmap[i])
    scatter!(ptemp_d, Js, χdwave, c=cmap[i], label=nothing)

    title!(ptemp_d, "d-wave")
    xlabel!(ptemp_d, "J")
    ylabel!(ptemp_d, "χ")
    title!(ptemp_s, "s-wave")
    xlabel!(ptemp_s, "J")
    ylabel!(ptemp_s, "χ")
end

# make a colourbar 
heatmap!(ptemp_s, zeros(2, 2), clims=(minimum(Ts), maximum(Ts)), cmap=:viridis, alpha=0)

if ndims == 3
    size_str = "$L × $L × $L"
elseif ndims == 2
    size_str = "$L × $L"
end

p2 = plot(ptemp_d, ptemp_s, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility for $size_str lattice with Q=$(round(Q, digits=3)), θ=$(θ_to_π(θ))")
if savefigs
    savefig(p2, joinpath(figpath, "susceptibility_temp_$(L)L.pdf"))
end
