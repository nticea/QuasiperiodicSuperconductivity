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

L = 11 # the full system is L × L 
ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
T_cutoff = 1e-1

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
df = load_dfs()
get_IPRs!(df) # in place!

# extract only the parameters we are interested in 
dfL = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q).&(df.ndims.==ndims), :]
dfL = dfL[(df.T.>=T_cutoff), :]

# use split-apply-combine
gdf = groupby(dfL, [:J, :T])
dfL = combine(gdf, [:χswave => mean, :χdwave => mean, :dχswave => mean, :dχdwave => mean,
    :IPR_swave_real => mean, :IPR_swave_momentum => mean, :IPR_x_real => mean, :IPR_x_momentum => mean,
    :IPR_y_real => mean, :IPR_y_momentum => mean, :IPR_z_real => mean, :IPR_z_momentum => mean])

px, pos = plot(margin=10Plots.mm), plot(margin=10Plots.mm)
pxd, posd = plot(margin=10Plots.mm), plot(margin=10Plots.mm)
Js = sort(unique(dfL.J))
cmap = cgrad(:matter, length(Js), categorical=true)
for (j, J) in enumerate(Js)
    dfJ = dfL[(dfL.J.==J), :]
    Ts = dfJ.T
    χswave = dfJ.χswave_mean
    χdwave = dfJ.χdwave_mean

    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χswave = χswave[sortidx]
    χdwave = χdwave[sortidx]

    # consider only the linear regime 
    # χs = χs[Ts.>=T_cutoff]
    # Ts = Ts[Ts.>=T_cutoff]

    if length(Ts) > 0
        plot!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(px, Ts, χdwave, color=cmap[j], label="J=$J, L=$L", xaxis=:log10)

        plot!(pos, Ts, χswave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(pos, Ts, χswave, color=cmap[j], label="J=$J, L=$L", xaxis=:log10)

        title!(px, "d-wave")
        xlabel!(px, "T")
        ylabel!(px, "χ")
        title!(pos, "s-wave")
        xlabel!(pos, "T")
        ylabel!(pos, "χ")
    end
end

ylims!(pos, (0, 7))
ylims!(px, (0, 10))
if ndims == 3
    size_str = "$L × $L × $L"
elseif ndims == 2
    size_str = "$L × $L"
end
p1 = plot(px, pos, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility for $size_str lattice with Q=$(round(Q, digits=3)), θ=$(θ_to_π(θ))")
if savefigs
    savefig(p1, joinpath(figpath, "susceptibility_J_$(L)L.pdf"))
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
    χswave = dfT.χswave_mean
    χdwave = dfT.χdwave_mean
    sortidx = sortperm(Js)
    Js = Js[sortidx]
    χswave = χswave[sortidx]
    χdwave = χdwave[sortidx]

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

p2 = plot(ptemp_d, ptemp_s, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Eigenvalue χ for $size_str lattice with Q=$(round(Q, digits=3)), θ=$(θ_to_π(θ))")
if savefigs
    savefig(p2, joinpath(figpath, "susceptibility_temp_$(L)L.pdf"))
end

px, pos = plot(margin=10Plots.mm), plot(margin=10Plots.mm)
pxd, posd = plot(margin=10Plots.mm), plot(margin=10Plots.mm)
Js = sort(unique(dfL.J))
cmap = cgrad(:matter, length(Js), categorical=true)
for (j, J) in enumerate(Js)
    dfJ = dfL[(dfL.J.==J), :]
    Ts = dfJ.T
    χswave = dfJ.dχswave_mean
    χdwave = dfJ.dχdwave_mean

    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χswave = χswave[sortidx]
    χdwave = χdwave[sortidx]

    # consider only the linear regime 
    # χs = χs[Ts.>=T_cutoff]
    # Ts = Ts[Ts.>=T_cutoff]

    if length(Ts) > 0
        plot!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(px, Ts, χdwave, color=cmap[j], label="J=$J, L=$L", xaxis=:log10)

        plot!(pos, Ts, χswave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(pos, Ts, χswave, color=cmap[j], label="J=$J, L=$L", xaxis=:log10)

        title!(px, "d-wave")
        xlabel!(px, "T")
        ylabel!(px, "χ")
        title!(pos, "s-wave")
        xlabel!(pos, "T")
        ylabel!(pos, "χ")
    end
end
ylims!(px, (-1, 0))
ylims!(pos, (-1, 0))


if ndims == 3
    size_str = "$L × $L × $L"
elseif ndims == 2
    size_str = "$L × $L"
end
p3 = plot(px, pos, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility for $size_str lattice with Q=$(round(Q, digits=3)), θ=$(θ_to_π(θ))")
if savefigs
    savefig(p1, joinpath(figpath, "susceptibility_J_$(L)L.pdf"))
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
    χswave = dfT.dχswave_mean
    χdwave = dfT.dχdwave_mean
    sortidx = sortperm(Js)
    Js = Js[sortidx]
    χswave = χswave[sortidx]
    χdwave = χdwave[sortidx]

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
#heatmap!(ptemp_s, zeros(2, 2), clims=(minimum(Ts), maximum(Ts)), cmap=:viridis, alpha=0)

p4 = plot(ptemp_d, ptemp_s, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Eigenvalue χ for $size_str lattice with Q=$(round(Q, digits=3)), θ=$(θ_to_π(θ))")
if savefigs
    savefig(p2, joinpath(figpath, "susceptibility_temp_$(L)L.pdf"))
end


## IPRS 

# also, plot the value of each J at various temperatures 
Ts = sort(unique(dfL.T))
T_min = Ts[1]
cmapreal = cgrad(:reds, 4, categorical=true)
cmapmomentum = cgrad(:blues, 4, categorical=true)
p_d = plot(margin=10Plots.mm, ylims=(-0.01, 1.1))
p_s = plot(margin=10Plots.mm, ylims=(-0.01, 1.1))
dfT = dfL[(dfL.T.==T_min), :]
Js = dfT.J

realIPR = dfT.IPR_swave_real_mean
kIPR = dfT.IPR_swave_momentum_mean
xrealIPR = dfT.IPR_x_real_mean
xmomentumIPR = dfT.IPR_x_momentum_mean
yrealIPR = dfT.IPR_x_real_mean
ymomentumIPR = dfT.IPR_x_momentum_mean
zrealIPR = dfT.IPR_x_real_mean
zmomentumIPR = dfT.IPR_x_momentum_mean

sortidx = sortperm(Js)
Js = Js[sortidx]
realIPR = realIPR[sortidx]
kIPR = kIPR[sortidx]
xrealIPR = xrealIPR[sortidx]
xmomentumIPR = xmomentumIPR[sortidx]
yrealIPR = yrealIPR[sortidx]
ymomentumIPR = ymomentumIPR[sortidx]
zrealIPR = zrealIPR[sortidx]
zmomentumIPR = zmomentumIPR[sortidx]

plot!(p_s, Js, realIPR, label=nothing, c=cmapreal[1])
scatter!(p_s, Js, realIPR, c=cmapreal[1], label="On-site real space")
plot!(p_d, Js, xrealIPR, label=nothing, c=cmapreal[2])
scatter!(p_d, Js, xrealIPR, c=cmapreal[2], label="x̂ real space")
plot!(p_d, Js, yrealIPR, label=nothing, c=cmapreal[3])
scatter!(p_d, Js, yrealIPR, c=cmapreal[3], label="ŷ real space")
plot!(p_d, Js, zrealIPR, label=nothing, c=cmapreal[4])
scatter!(p_d, Js, zrealIPR, c=cmapreal[4], label="ẑ real space")

plot!(p_s, Js, kIPR, label=nothing, c=cmapmomentum[1])
scatter!(p_s, Js, kIPR, c=cmapmomentum[1], label="On-site momentum space")
plot!(p_d, Js, xmomentumIPR, label=nothing, c=cmapmomentum[2])
scatter!(p_d, Js, xmomentumIPR, c=cmapmomentum[2], label="x̂ momentum space")
plot!(p_d, Js, ymomentumIPR, label=nothing, c=cmapmomentum[3])
scatter!(p_d, Js, ymomentumIPR, c=cmapmomentum[3], label="ŷ momentum space")
plot!(p_d, Js, zmomentumIPR, label=nothing, c=cmapmomentum[4])
scatter!(p_d, Js, zmomentumIPR, c=cmapmomentum[4], label="ẑ momentum space")

title!(p_d, "d-wave")
xlabel!(p_d, "J")
ylabel!(p_d, "IPR")
title!(p_s, "s-wave")
xlabel!(p_s, "J")
ylabel!(p_s, "IPR")


# # make a colourbar 
# heatmap!(p_real, zeros(1, 1), clims=(minimum(Ts), maximum(Ts)), cmap=:viridis, alpha=0)

p5 = plot(p_s, p_d, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="IPR for $size_str lattice with Q=$(round(Q, digits=3)), θ=$(θ_to_π(θ))")
