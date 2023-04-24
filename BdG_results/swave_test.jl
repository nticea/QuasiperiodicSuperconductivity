## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
include("../src/BdG.jl")

using Plots
using ProgressBars
using CSV
using DataFrames

## PARAMETERS ##
L = 55 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
θ = nothing #π / 7
μ = 0
pairing_symmetry = "s-wave"
tol = 1e-4

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_results.csv")

# J, V0, T 
Js = [0, 1, 2, 3]
V0s = expspace(-0.8, 0.7, 20)
Ts = expspace(-3, 0, 20)
Δs = zeros(length(Js), length(V0s), length(Ts))

# load in the dataframe, if it exists. If not, make a new one
df = load_dataframe(savepath)

# the finite size gap is determined by the largest gap in the J=0 system 
ΔE = finite_size_gap(L=L, t=t, J=0, Q=Q, μ=μ, θ=θ)
fsgap = maximum(ΔE)
@show fsgap

## RUNNING THE CODE ## 
for J in Js # iterate through all J values
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))
    for i in iter # iterate through all temperatures
        T = Ts[i]
        for V0 in V0s # iterature through all V0 values 
            # check whether this particular (J,T,V0) combo has been already computed 
            if !already_calculated(df; L=L, J=J, V0=V0, T=T)
                # calculate Δ at a given (J,T,V0)
                Δ = compute_Δ(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol, θ=θ)
                update_results!(df; L=L, λ=Δ, J=J, V0=V0, T=T)
                CSV.write(savepath, df)
                flush(stdout)
            end
        end
    end
end

# # plot some things 
# plot(conv)
# title!("Convergence for J=$J (L=$L)")
# xlabel!("Iteration number")

# Δmat = reshape(Δi, L, L)
# heatmap(Δmat)
# xticks!(collect(1:5:L))
# yticks!(collect(1:5:L))
# xlabel!("Site (x)")
# ylabel!("Site (y)")
# title!("Gap (J=$J, T=$T, V0=$V0)")
