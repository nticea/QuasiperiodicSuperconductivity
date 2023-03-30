## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using ProgressBars
using CSV
using DataFrames

## PARAMETERS ##
L = 49 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
pairing_symmetry = "s-wave"

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_results_J1.csv")

# J, V0, T 
Js = [1]
V0s = expspace(-0.5, 1.5, 20)
Ts = expspace(-3.5, 0, 20)
λs = zeros(length(Js), length(V0s), length(Ts))

# load in the dataframe, if it exists. If not, make a new one
df = load_dataframe(savepath)

## RUNNING THE CODE ## 
for J in Js # iterate through all J values
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))
    for t in iter # iterate through all temperatures
        T = Ts[t]
        for V0 in V0s # iterature through all V0 values 
            # check whether this particular (J,T,V0) combo has been already computed 
            if !already_calculated(df; L=L, J=J, V0=V0, T=T)
                # calculate λmax at a given (J,T,V0)
                λ = λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
                update_results!(df; L=L, λ=λ, J=J, V0=V0, T=T)
                CSV.write(savepath, df)
                flush(stdout)
            end
        end
    end
end