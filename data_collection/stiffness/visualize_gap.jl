## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using CSV
using DataFrames
using StatsPlots, Statistics

include("../../src/BdG.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")
include("utilities.jl")

## PARAMETERS ## 

L = 7 # the full system is L × L 
t = 1
Q = (√5 - 1) / 2
μ = 0.75
θ = π / 7
V0 = -3#2#-3
V1 = 0#-1.5#0
ndims = 3
periodic = true
disorder = false
savefigs = false
slice = 3

if disorder
    dirname = "data_$(ndims)D_disorder"
else
    dirname = "data_$(ndims)D_QP"
end
if ndims == 3
    size_str = "$L × $L × $L"
elseif ndims == 2
    size_str = "$L × $L"
end

# read files 
if savefigs
    mkpath(joinpath(@__DIR__, "figures"))
end
# df_LGE_full = load_LGE(dirname)
# df_BdG_full = load_BdG(dirname)

# extract only the parameters we are interested in 
df_LGE = df_LGE_full[(df_LGE_full.L.==L).&(df_LGE_full.ndims.==ndims).&(df_LGE_full.θ.==θ).&(df_LGE_full.Q.==Q).&(df_LGE_full.V0.==V0).&(df_LGE_full.V1.==V1).&(df_LGE_full.μ.==μ), :]
df_BdG = df_BdG_full[(df_BdG_full.L.==L).&(df_BdG_full.ndims.==ndims).&(df_BdG_full.θ.==θ).&(df_BdG_full.Q.==Q).&(df_BdG_full.V0.==V0).&(df_BdG_full.V1.==V1).&(df_BdG_full.μ.==μ), :]

df_BdG0 = df_BdG[(df_BdG.T.==0), :]
gdf = groupby(df_BdG0, [:J])
g = gdf[1]
Δsavg = zeros(size(g.Δ[1])) .* 1im
for Δ in g.Δ
    Δsavg .+= Δ
end
@assert 1 == 0

J = 1
dfsub = df_LGE[(df_LGE.J.==J), :]
ϕx, ϕy, ϕz = dfsub.ϕx[1], dfsub.ϕy[1], dfsub.ϕz[1]
dfsub = dfsub[(dfsub.J.==J).&(dfsub.ϕx.==ϕx).&(dfsub.ϕy.==ϕy).&(dfsub.ϕz.==ϕz).&(abs.(dfsub.λ .- 1).<1e-2), :]
m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
Δ_LGE = dfsub.Δ[1]
Tc = dfsub.T[1]
Δ = spatial_profile(m, Δ=abs.(Δ_LGE))
evs = abs.(Δ)[:, :, :, slice]

p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false, aspect_ratio=true)
for x in 1:L
    for y in 1:L

        # bonds 
        plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)
        plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[5, x, y]), alpha=10 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)

        # onsite dot 
        scatter!(p, [x], [y], ms=100 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
    end
end

xlabel!(p, "Site (x)")
ylabel!(p, "Site, (y)")
title!(p, "J=$(m.J), V0=$(m.V0), V1=$(m.V1), θ=$(θ_to_π(m.θ)) on $(m.L)×$(m.L)×$(m.L) lattice at z=$slice")

# J = 1
# dfsub = df_BdG[(df_BdG.J.==J), :]
# ϕx, ϕy, ϕz = dfsub.ϕx[1], dfsub.ϕy[1], dfsub.ϕz[1]
# dfsub = dfsub[(dfsub.J.==J).&(dfsub.ϕx.==ϕx).&(dfsub.ϕy.==ϕy).&(dfsub.ϕz.==ϕz).&(dfsub.T.==0), :]
# m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
# Δ_LGE = dfsub.Δ[1]
# Tc = dfsub.T[1]
# Δ = spatial_profile(m, Δ=abs.(Δ_LGE))
# evs = abs.(Δ)[:, :, :, slice]

# p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false, aspect_ratio=true)
# for x in 1:L
#     for y in 1:L

#         # bonds 
#         plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
#         plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
#         plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)
#         plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[5, x, y]), alpha=10 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)

#         # onsite dot 
#         scatter!(p, [x], [y], ms=100 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
#     end
# end

# xlabel!(p, "Site (x)")
# ylabel!(p, "Site, (y)")
# title!(p, "J=$(m.J), V0=$(m.V0), V1=$(m.V1), θ=$(θ_to_π(m.θ)) on $(m.L)×$(m.L)×$(m.L) lattice at z=$slice")


@assert 1 == 0

# Spatial profiles of LGE soln at Tc (real space and configuration space)
ϕx, ϕy, ϕz = df_LGE.ϕx[1], df_LGE.ϕy[1], df_LGE.ϕz[1]
Js_to_plot = [0.25]
ps_config = []
for J in Js_to_plot
    dfsub = df_LGE[(df_LGE.J.==J).&(df_LGE.ϕx.==ϕx).&(df_LGE.ϕy.==ϕy).&(df_LGE.ϕz.==ϕz).&(abs.(df_LGE.λ .- 1).<1e-2), :]
    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    if size(dfsub)[1] > 0
        if V1 != 0
            Δ_LGE = dfsub.Δ[1]
            #χ = symmetry_character(m, Δ=Δ_LGE)
            Tc = dfsub.T[1]
            Δ = spatial_profile(m, Δ=real.(Δ_LGE))
            p0 = plot_spatial_profile(m, Δ=real.(Δ_LGE), title="J=$J at Tc=$(round(Tc,digits=2))")
            # take a slice if in 3D 
            if ndims == 3
                Δ = Δ[:, :, :, slice]
            end
            ps = plot_in_configuration_space(m, Δ=abs.(Δ))
            push!(ps_config, ps)
            if savefigs
                savefig(p, joinpath(@__DIR__, "figures", "spatial_profile_$(L)L_$(V0)V0_$(V1)V1_J$J.pdf"))
            end
        else
            Δ_LGE = dfsub.Δ[1]
            Tc = dfsub.T[1]
            Δ = spatial_profile(m, Δ=Δ_LGE)
            # take a slice if in 3D 
            if ndims == 3
                Δ = Δ[:, :, :, slice]
            end
            p = plot_in_configuration_space(m, Δ=Δ)
            push!(ps_config, p)
            if savefigs
                savefig(p, joinpath(@__DIR__, "figures", "spatial_profile_$(L)L_$(V0)V0_$(V1)V1_J$J.pdf"))
            end
        end
    end
end

# Spatial profiles of LGE soln at Tc (real space only)
# Js_to_plot = [0, 0.6, 0.7, 0.9, 1, 1.1, 1.2, 1.5, 2, 2.5, 3]
# ps_real = []
# for J in Js_to_plot
#     dfsub = df_LGE_Tc[(df_LGE_Tc.L.==L).&(df_LGE_Tc.J.==J).&(df_LGE_Tc.θ.==θ).&(df_LGE_Tc.ϕx.==ϕx).&(df_LGE_Tc.ϕy.==ϕy).&(df_LGE_Tc.Q.==Q).&(df_LGE_Tc.V0.==V0).&(df_LGE_Tc.V1.==V1), :]
#     if size(dfsub)[1] > 0
#         m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
#         Δ_LGE = dfsub.Δ[1]
#         Tc = dfsub.T[1]
#         p0 = plot_spatial_profile(m, Δ=Δ_LGE, title="J=$J at Tc=$(round(Tc,digits=2))")
#         push!(ps_real, p0)
#     end
# end

# p = plot(ps_real..., layout=Plots.grid(2, 6,
#         widths=[1 / 6 for _ in 1:6]), size=(1500, 800), aspect_ratio=:equal, plot_title=" LGE Δ(V0=$V0, V1=$V1, θ=$(θ_to_π(θ)), ϕx=$(θ_to_π(ϕx)), ϕy=$(θ_to_π(ϕy))) for $size_str lattice")

# if savefigs
#     savefig(p, joinpath(@__DIR__, "figures", "real_spatial_profile_$(L)L_$(V0)V0_$(V1)V1.pdf"))
# end