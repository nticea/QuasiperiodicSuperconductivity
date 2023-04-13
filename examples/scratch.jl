## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
include("../src/meanfield.jl")

using Plots
using ProgressBars
using CSV
using DataFrames

## PARAMETERS ##
L = 9 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 1e-8
pairing_symmetry = "s-wave"

# saving information 
savepath = joinpath(@__DIR__, "$(L)Nx$(L)Ny_results.csv")

# J, V0, T 
Js = [0, 1, 2, 3]
V0s = expspace(-0.8, 0.7, 20)
Ts = expspace(-3, 0, 20)
λs = zeros(length(Js), length(V0s), length(Ts))

# load in the dataframe, if it exists. If not, make a new one
df = load_dataframe(savepath)

## RUNNING THE CODE ## 
for (k, J) in enumerate(Js)
    println("Running J=$(J)")
    iter = ProgressBar(1:length(Ts))
    for i in iter # iterate through all temperatures
        T = Ts[i]
        for (j, V0) in enumerate(V0s) # iterature through all V0 values 
            # check whether this particular (J,T,V0) combo has been already computed 
            if !already_calculated(df; L=L, J=J, V0=V0, T=T)
                # calculate λmax at a given (J,T,V0)
                λ = λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
                update_results!(df; L=L, λ=λ, J=J, V0=V0, T=T)
                CSV.write(savepath, df)
            end
        end
    end
end

λs = zeros(length(Js), length(V0s), length(Ts))
for (k, J) in enumerate(Js)
    for (i, T) in enumerate(Ts)
        for (j, V0) in enumerate(V0s)
            λ = λmax(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0)
            @show J, T, V0, λ
            λs[k, j, i] = λ
        end
    end
end

# do interpolations 
Tcs = zeros(length(Js), length(V0s))
for k in 1:length(Js)
    for j in 1:length(V0s)
        knots = reverse(λs[k, j, :])
        Interpolations.deduplicate_knots!(knots)
        try
            interp_linear = linear_interpolation(knots, reverse(Ts))
            Tcs[k, j] = interp_linear(1)
        catch e
            Tcs[k, j] = NaN
        end
    end
end

p2 = plot()
cmap = cgrad(:Set1_9, length(V0s), categorical=true)
for (k, J) in enumerate(Js)
    plot!(p2, V0s, Tcs[k, :], xaxis=:log10, yaxis=:log10, c=cmap[k], label=nothing)
    scatter!(p2, V0s, Tcs[k, :], xaxis=:log10, yaxis=:log10, c=cmap[k], label="J=$(J)")
end
title!(p2, "Transition temperature for $(L)x$(L) square lattice")
xlabel!(p2, "V")
ylabel!(p2, "Tc")


