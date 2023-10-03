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

L = 17 # the full system is L × L 
ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
T_cutoff = 0#1e-1

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
dirname = "data_PBC"
df = load_dfs(dirname=dirname)

function mean_χ(χs)
    χs = hcat(χs...)
    return mean(χs, dims=2)
end

px, pos = plot(margin=10Plots.mm), plot(margin=10Plots.mm)
# extract only the parameters we are interested in 
df = df[(df.T.>=T_cutoff).&(df.L.==L), :]
df = df[(df.ϕx.==df.ϕx[1]).&(df.ϕy.==df.ϕy[1]).&(df.ϕz.==df.ϕz[1]), :]


gdf = groupby(df, [:J, :T])
for g in gdf
    # replace χ and dχlogT with their means 
    g.χ .= [mean_χ(g.χ) for _ in 1:length(g.χ)]
    g.dχdlogT .= [mean_χ(g.dχdlogT) for _ in 1:length(g.dχdlogT)]
end
df = vcat(gdf...)

# compute the eigenvector at a fixed temperature (T_cutoff) for each J 
gdf = groupby(df, [:J])
for g in gdf
    lg = size(g)[1]
    # find the temperature closest to T_cutoff 
    idx = closest_indices(g.T, T_cutoff, 1)[1]
    # get the eigenvector of the susceptibility
    χ = reshape(g.χ[idx], (4, 4))
    λs, Us = eigen(χ[2:end, 2:end], sortby=real)
    g.χ_λ = [λs[end] for _ in 1:lg]
    g.χ_U = [Us[:, end] for _ in 1:lg]
end
df = vcat(gdf...)
dfL = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q).&(df.ndims.==ndims), :]
Js = sort(unique(dfL.J))
cmap = cgrad(:matter, length(Js), categorical=true)

for (j, J) in enumerate(Js)
    dfJ = dfL[(dfL.J.==J), :]
    Ts = dfJ.T
    χs = dfJ.χ
    Us = dfJ.χ_U

    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χs = χs[sortidx]
    Us = Us[sortidx]
    χs = [reshape(χ, 4, 4) for χ in χs]

    # consider only the linear regime 
    # χs = χs[Ts.>=T_cutoff]
    # Ts = Ts[Ts.>=T_cutoff]

    if length(Ts) > 0

        χswave, χdwave = [], []
        for (χ, U) in zip(χs, Us)
            χsw, χdw = uniform_susceptibility_components(χ, U=U, ndims=ndims)
            push!(χswave, χsw)
            push!(χdwave, χdw)
        end

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

if ndims == 3
    size_str = "$L × $L × $L"
elseif ndims == 2
    size_str = "$L × $L"
end
p1 = plot(px, pos, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility for $size_str lattice with Q=$(round(Q, digits=3))
        ")
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
    χs = dfT.χ
    Us = dfT.χ_U

    sortidx = sortperm(Js)
    Js = Js[sortidx]
    χs = χs[sortidx]
    Us = Us[sortidx]
    χs = [reshape(χ, 4, 4) for χ in χs]

    χswave, χdwave = [], []
    for (χ, U) in zip(χs, Us)
        χsw, χdw = uniform_susceptibility_components(χ, U=U, ndims=ndims)
        push!(χswave, χsw)
        push!(χdwave, χdw)
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
#heatmap!(ptemp_s, zeros(2, 2), clims=(minimum(Ts), maximum(Ts)), cmap=:viridis, alpha=0)

p2 = plot(ptemp_d, ptemp_s, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="Susceptibility for $size_str lattice with Q=$(round(Q, digits=3))
        ")
if savefigs
    savefig(p2, joinpath(figpath, "susceptibility_temp_$(L)L.pdf"))
end

px, pos = plot(margin=10Plots.mm, ylims=(-1, 1)), plot(margin=10Plots.mm, ylims=(-1, 1))
# extract only the parameters we are interested in 
dfL = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q).&(df.ndims.==ndims), :]
Js = sort(unique(dfL.J))
cmap = cgrad(:matter, length(Js), categorical=true)

# normalization factor 
ndf = dfL[(dfL.J.==0).&&(dfL.T.==minimum(dfL.T)), :]
if size(ndf)[1] > 0
    nfχ = ndf.dχdlogT[1]
    U = ndf.χ_U[1]
    nfs, nfd = 1, 1#uniform_susceptibility_components(nfχ, U=U, ndims=ndims)
else
    nfs, nfd = 1, 1
end

for (j, J) in enumerate(Js)
    dfJ = dfL[(dfL.J.==J), :]
    Ts = dfJ.T
    χs = dfJ.dχdlogT
    Us = dfJ.χ_U

    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χs = χs[sortidx]
    Us = Us[sortidx]
    χs = [reshape(χ, 4, 4) for χ in χs]

    # consider only the linear regime 
    # χs = χs[Ts.>=T_cutoff]
    # Ts = Ts[Ts.>=T_cutoff]

    if length(Ts) > 0
        χswave, χdwave = [], []
        for (χ, U) in zip(χs, Us)
            χsw, χdw = uniform_susceptibility_components(χ, U=U, ndims=ndims)
            push!(χswave, χsw / abs(nfs))
            push!(χdwave, χdw / abs(nfd))
        end

        plot!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(px, Ts, χdwave, color=cmap[j], label=nothing, xaxis=:log10)

        plot!(pos, Ts, χswave, color=cmap[j], label=nothing, xaxis=:log10)
        scatter!(pos, Ts, χswave, color=cmap[j], label="J=$J", xaxis=:log10)

        title!(px, "d-wave")
        xlabel!(px, "T")
        ylabel!(px, "dχdlogT")
        title!(pos, "s-wave")
        xlabel!(pos, "T")
        ylabel!(pos, "dχdlogT")
    end
end
ylims!(pos, (-0.3, 0))
ylims!(px, (-0.3, 0))

if ndims == 3
    size_str = "$L × $L × $L"
elseif ndims == 2
    size_str = "$L × $L"
end
p3 = plot(px, pos, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="dχdlogT for $size_str lattice with Q=$(round(Q, digits=3))
        ")
if savefigs
    savefig(p3, joinpath(figpath, "dχdlogT_$(L)L.pdf"))
end

Ts = sort(unique(dfL.T))
Ts = Ts[Ts.>=T_cutoff]
cmap = cgrad(:viridis, length(Ts), categorical=true)

# finally, plot the dχdlogT data
ptemp_s = plot(margin=10Plots.mm)
ptemp_d = plot(margin=10Plots.mm)

for (i, T) in enumerate(Ts)
    dfT = dfL[(dfL.T.==T), :]
    Js = dfT.J
    χs = dfT.dχdlogT
    Us = dfT.χ_U

    sortidx = sortperm(Js)
    Js = Js[sortidx]
    χs = χs[sortidx]
    Us = Us[sortidx]
    χs = [reshape(χ, 4, 4) for χ in χs]

    χswave, χdwave = [], []
    for (χ, U) in zip(χs, Us)
        χsw, χdw = uniform_susceptibility_components(χ, U=U, ndims=ndims)
        push!(χswave, χsw / abs(nfs))
        push!(χdwave, χdw / abs(nfd))
    end

    plot!(ptemp_s, Js, χswave, label=nothing, c=cmap[i])
    scatter!(ptemp_s, Js, χswave, c=cmap[i], label=nothing)

    plot!(ptemp_d, Js, χdwave, label=nothing, c=cmap[i])
    scatter!(ptemp_d, Js, χdwave, c=cmap[i], label=nothing)

    title!(ptemp_d, "d-wave")
    xlabel!(ptemp_d, "J")
    ylabel!(ptemp_d, "dχ/dlog₁₀T")
    title!(ptemp_s, "s-wave")
    xlabel!(ptemp_s, "J")
    ylabel!(ptemp_s, "dχ/dlog₁₀T")
end

# make a colourbar 
#heatmap!(ptemp_s, zeros(2, 2), clims=(minimum(Ts), maximum(Ts)), cmap=:viridis, alpha=0)

p4 = plot(ptemp_d, ptemp_s, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="dχ/dlog₁₀T for $size_str lattice with Q=$(round(Q, digits=3))")
if savefigs
    savefig(p4, joinpath(figpath, "susceptibility_dχdlogT_$(L)L.pdf"))
end

# plot the symmetry 
gdf = groupby(dfL, [:J])
Js = []
x1, x2, x3 = [], [], []
λs = []
for g in gdf
    push!(Js, g.J[1])
    U = g.χ_U[1]
    if real(U[1]) < 0
        U *= -1
    end
    push!(x1, U[1])
    push!(x2, U[2])
    push!(x3, U[3])
    push!(λs, g.χ_λ[1])
end

p5 = plot(Js, real.(x1), color="red", label=nothing)
scatter!(p5, Js, real.(x1), color="red", label="real")
plot!(p5, Js, real.(x2), color="blue", label=nothing)
scatter!(p5, Js, real.(x2), color="blue", label="real")
plot!(p5, Js, real.(x3), color="green", label=nothing)
scatter!(p5, Js, real.(x3), color="green", label="real")
plot!(Js, imag.(x1), color="red", label="imaginary", ls=:dashdot)
plot!(p5, Js, imag.(x2), color="blue", label="imaginary", ls=:dashdot)
plot!(p5, Js, imag.(x3), color="green", label="imaginary", ls=:dashdot)
xlabel!(p5, "J")
ylabel!(p5, "Loading onto component")

p6 = plot(Js, real.(λs), color="red", label=nothing)
scatter!(p6, Js, real.(λs), color="red", label="real part")
plot!(p6, Js, imag.(λs), color="blue", label=nothing)
scatter!(p6, Js, imag.(λs), color="blue", label="imaginary part")
xlabel!(p6, "J")
ylabel!(p6, "λ")
title!("λ(χ)")