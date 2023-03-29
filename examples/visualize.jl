## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
include("../src/model.jl")
using Plots

savepath = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/29Nx29Ny_results_J1.h5"
savepath_df = "/Users/nicole/Dropbox/Grad school/Trithep/quasiperiodic/QuasiperiodicSuperconductivity/examples/29Nx29Ny_results.csv"
results = load_results(savepath)
plot_Tcs(results)

L, Js, Ts, V0s, λs = results.L, results.Js, results.Ts, results.V0s, results.λs

Jarr = zeros(length(Js), length(V0s), length(Ts))
Tarr = zeros(length(Js), length(V0s), length(Ts))
Varr = zeros(length(Js), length(V0s), length(Ts))
Larr = L * ones(length(Js), length(V0s), length(Ts))

for (k, J) in enumerate(Js) # iterate through all J values
    for (i, T) in enumerate(Ts) # iterate through all temperatures
        for (j, V0) in enumerate(V0s) # iterature through all V0 values 
            # calculate λmax at a given (J,T,V0)
            Jarr[k, j, i] = J
            Tarr[k, j, i] = T
            Varr[k, j, i] = V0
        end
    end
end

# reshape everything 
Jarr = reshape(Jarr, prod(size(Jarr)))
Tarr = reshape(Tarr, prod(size(Tarr)))
Varr = reshape(Varr, prod(size(Varr)))
λs = reshape(λs, prod(size(λs)))
Ls = L * ones(length(λs))

df1 = DataFrame(L=Ls, T=Tarr, V0=Varr, λ=λs, J=Jarr)
#CSV.write(savepath_df, df1)