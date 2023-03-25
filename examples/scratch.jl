## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots
using ProgressBars

## PARAMETERS ##
L = 29 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
pairing_symmetry = "s-wave"

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_results.npy")

# J, V0, T 
Js = [0, 1, 2, 3]
V0s = expspace(-0.8, 0.7, 20)
Ts = expspace(-3, 0, 20)
λs = zeros(length(Js), length(V0s), length(Ts))

# for storing the results
global results = Results(L, λs, Js, V0s, Ts)

## RUNNING THE CODE ## 
for (k, J) in enumerate(Js) # iterate through all J values
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))
    for i in iter # iterate through all temperatures
        T = Ts[i]
        for (j, V0) in enumerate(V0s) # iterature through all V0 values 
            # calculate λmax at a given (J,T,V0)
            λs[k, j, i] = λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
        end

        # interim update results & save
        global results = Results(L, λs, Js, V0s, Ts)
        save_structs(results, savepath)
    end
end

# interpolate to find the Tcs and then plot 
plot_Tcs(results)