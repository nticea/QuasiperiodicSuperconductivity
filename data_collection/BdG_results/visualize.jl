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
V0 = 1
V1 = -1.5

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

# Spatial profiles of LGE and BdG solns at Tc (real space and configuration space)
ϕx, ϕy = 0.123, 0.987
J = 2.5
dfsub = df_LGE_Tc[(df_LGE_Tc.L.==L).&(df_LGE_Tc.J.==J).&(df_LGE_Tc.θ.==θ).&(df_LGE_Tc.ϕx.==ϕx).&(df_LGE_Tc.ϕy.==ϕy).&(df_LGE_Tc.Q.==Q).&(df_LGE_Tc.V0.==V0).&(df_LGE_Tc.V1.==V1), :]
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
#savefig(p, joinpath(@__DIR__, "figures", "spatial_profile_J$J.pdf"))

h = plot_potential(L=L, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy)