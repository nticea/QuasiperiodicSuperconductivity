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

L = 11 # the full system is L × L 
Q = (√5 - 1) / 2
θ = π / 7
V0 = -2.3
V1 = 0
savefigs = false

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
df_LGE = df_LGE[(df_LGE.L.==L).&(df_LGE.θ.==θ).&(df_LGE.Q.==Q).&(df_LGE.V0.==V0).&(df_LGE.V1.==V1), :]
df_BdG = df_BdG[(df_BdG.L.==L).&(df_BdG.θ.==θ).&(df_BdG.Q.==Q).&(df_BdG.V0.==V0).&(df_BdG.V1.==V1), :]

## FIGURES ##
# Comparing Tc from LGE and from BdG 
df_LGE_Tc = df_LGE[df_LGE.T.>0, :]
df_BdG_Tc = df_BdG[df_BdG.T.>0, :]

p1 = plot()
grouped = groupby(df_LGE_Tc, [:ϕx, :ϕy])
cmap = cgrad(:matter, length(grouped), categorical=true)
for (i, df) in enumerate(grouped)
    ϕx, ϕy = unique(df.ϕx)[1], unique(df.ϕy)[1]
    Js, Tcs = df.J, df.T
    sortidx = sortperm(Js)
    Js, Tcs = Js[sortidx], Tcs[sortidx]
    plot!(p1, Js, Tcs, color=cmap[i], label=nothing)
    scatter!(p1, Js, Tcs, color=cmap[i], label="(ϕx,ϕy)=($(round(ϕx, digits=2)),$(round(ϕy, digits=2)))", ylabel="Tc",)
end
xlabel!(p1, "J")
title!(p1, "Tc for V0=$V0, V1=$V1, θ=$(θ_to_π(θ))\n on $L × $L lattice")

# Superfluid stiffness vs Tc 
df_LGE_0 = df_LGE[df_LGE.T.==0, :]
df_BdG_0 = df_BdG[df_BdG.T.==0, :]

Js = sort(unique(df_BdG_0.J))
Ds = zeros(length(Js), 4)
for (j, J) in enumerate(Js)
    dfsub = df_BdG_0[(df_BdG_0.L.==L).&(df_BdG_0.J.==J).&(df_BdG_0.θ.==θ).&(df_BdG_0.Q.==Q).&(df_BdG_0.V0.==V0).&(df_BdG_0.V1.==V1), :]
    Ks, Πs = hcat(dfsub.K...), hcat(dfsub.Π...)
    Ks, Πs = mean(Ks, dims=2), mean(Πs, dims=2)
    if J == 0
        @show V0, V1
        @show Ks, Πs
    end
    Ds[j, :] = -Ks + Πs
end

dirs = ["Dₛ/π (x̂)", "Dₛ/π (ŷ)"]
cmap = ["green", "orange"]
for i in 1:2
    plot!(p1, Js, Ds[:, i], label=nothing, c=cmap[i], secondary=true)
    scatter!(p1, Js, Ds[:, i], label=dirs[i], c=cmap[i], secondary=true)
end
plot!(p1, legend=:right)

if savefigs
    savefig(p1, joinpath(@__DIR__, "figures", "$(L)L_$(V0)V0_$(V1)V1_stiffness_averaged.pdf"))
end

p2 = plot()
grouped = groupby(df_BdG_0, [:ϕx, :ϕy])
cmap = cgrad(:matter, length(grouped), categorical=true)
for (i, df) in enumerate(grouped)
    ϕx, ϕy = unique(df.ϕx)[1], unique(df.ϕy)[1]
    Js, Ks, Πs = df.J, df.K, df.Π
    sortidx = sortperm(Js)
    Js, Ks, Πs = Js[sortidx], Ks[sortidx], Πs[sortidx]
    Ds = -Ks + Πs
    Ds = [D[2] for D in Ds]
    plot!(p2, Js, Ds, color=cmap[i], label=nothing)
    scatter!(p2, Js, Ds, color=cmap[i], label="(ϕx,ϕy)=($(round(ϕx, digits=2)),$(round(ϕy, digits=2)))")
end
xlabel!(p2, "J")
ylabel!(p2, "Ds/π")
title!(p2, "Tc₂ for V0=$V0, V1=$V1, θ=$(θ_to_π(θ))\n on $L × $L lattice")
if savefigs
    savefig(p2, joinpath(@__DIR__, "figures", "$(L)L_$(V0)V0_$(V1)V1_stiffness_unaveraged.pdf"))
end

# Spatial profiles of LGE soln at Tc (real space and configuration space)
ϕx, ϕy = 0, 0
Js_to_plot = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 2.5, 3, 3.5]
ps_config = []
for J in Js_to_plot
    dfsub = df_LGE_Tc[(df_LGE_Tc.L.==L).&(df_LGE_Tc.J.==J).&(df_LGE_Tc.θ.==θ).&(df_LGE_Tc.ϕx.==ϕx).&(df_LGE_Tc.ϕy.==ϕy).&(df_LGE_Tc.Q.==Q).&(df_LGE_Tc.V0.==V0).&(df_LGE_Tc.V1.==V1), :]
    if size(dfsub)[1] == 1
        if V1 != 0
            Δ_LGE = dfsub.Δ[1]
            χ = symmetry_character(Δ_LGE, L=L)
            Tc = dfsub.T[1]
            Δ = spatial_profile(Δ_LGE, L=L)
            p0 = plot_spatial_profile(Δ, L=L, title="Real space, symmetry=$(round(χ,digits=2))")
            p1 = plot_in_config_space(Δ[5, :, :], L=L, Q=Q, θ=θ, title="On-site")
            p2 = plot_in_config_space(Δ[1, :, :], L=L, Q=Q, θ=θ, title="-x̂ bond")
            p3 = plot_in_config_space(Δ[2, :, :], L=L, Q=Q, θ=θ, title="+ŷ bond")
            p4 = plot_in_config_space(Δ[3, :, :], L=L, Q=Q, θ=θ, title="x̂ bond")
            p5 = plot_in_config_space(Δ[4, :, :], L=L, Q=Q, θ=θ, title="-ŷ bond")
            p = plot(p1, p2, p3, p4, p5, p0, layout=Plots.grid(2, 3,
                    widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1000), aspect_ratio=:equal, plot_title=" LGE Δ(J=$J, V0=$V0, V1=$V1, θ=$(θ_to_π(θ)), ϕx=$(θ_to_π(ϕx)), ϕy=$(θ_to_π(ϕy))) for $L × $L lattice at Tc=$(round(Tc,digits=2))")
            push!(ps_config, p)
            if savefigs
                savefig(p, joinpath(@__DIR__, "figures", "spatial_profile_$(L)L_$(V0)V0_$(V1)V1_J$J.pdf"))
            end
        else
            Δ_LGE = dfsub.Δ[1]
            Tc = dfsub.T[1]
            Δ = spatial_profile(Δ_LGE, L=L)
            χ = symmetry_character(Δ_LGE, L=L)
            p0 = plot_spatial_profile(Δ, L=L, title="Real space, symmetry=$(round(χ,digits=2))")
            p1 = plot_in_config_space(Δ[5, :, :], L=L, Q=Q, θ=θ, title="On-site")
            p = plot(p1, p0, layout=Plots.grid(1, 2,
                    widths=[1 / 2, 1 / 2]), size=(1000, 700), aspect_ratio=:equal, plot_title=" LGE Δ(J=$J, V0=$V0, V1=$V1, θ=$(θ_to_π(θ)), ϕx=$(θ_to_π(ϕx)), ϕy=$(θ_to_π(ϕy))) for $L × $L lattice at Tc=$(round(Tc,digits=2))")
            push!(ps_config, p)
            if savefigs
                savefig(p, joinpath(@__DIR__, "figures", "spatial_profile_$(L)L_$(V0)V0_$(V1)V1_J$J.pdf"))
            end
        end

    end
end

# Spatial profiles of LGE soln at Tc (real space only)
ϕx, ϕy = 0, 0
Js_to_plot = [0, 0.5, 0.7, 0.9, 1, 1.1, 1.2, 1.5, 2, 2.5]
ps_real = []
for J in Js_to_plot
    dfsub = df_LGE_Tc[(df_LGE_Tc.L.==L).&(df_LGE_Tc.J.==J).&(df_LGE_Tc.θ.==θ).&(df_LGE_Tc.ϕx.==ϕx).&(df_LGE_Tc.ϕy.==ϕy).&(df_LGE_Tc.Q.==Q).&(df_LGE_Tc.V0.==V0).&(df_LGE_Tc.V1.==V1), :]
    if size(dfsub)[1] == 1
        Δ_LGE = dfsub.Δ[1]
        Tc = dfsub.T[1]
        Δ = spatial_profile(Δ_LGE, L=L)
        p0 = plot_spatial_profile(Δ, L=L, title="J=$J at Tc=$(round(Tc,digits=2))")
        push!(ps_real, p0)
    end
end

p = plot(ps_real..., layout=Plots.grid(2, 5,
        widths=[1 / 5 for _ in 1:5]), size=(1500, 800), aspect_ratio=:equal, plot_title=" LGE Δ(V0=$V0, V1=$V1, θ=$(θ_to_π(θ)), ϕx=$(θ_to_π(ϕx)), ϕy=$(θ_to_π(ϕy))) for $L × $L lattice")

if savefigs
    savefig(p, joinpath(@__DIR__, "figures", "real_spatial_profile_$(L)L_$(V0)V0_$(V1)V1.pdf"))
end