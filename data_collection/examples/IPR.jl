## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/model.jl")
include("../../src/results.jl")

## PARAMETERS ## 
L = 11
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 1.5
V1 = -1
ϕx = 0
ϕy = 0
ϕz = 0
periodic = true
ndims = 3
nbins = 41

Js = collect(0:0.25:3)#expspace(log10(2) - 1, log10(2) + 1, 30)
iprs_real = zeros(length(Js), nbins)
iprs_k = zeros(length(Js), nbins)

dos = zeros(length(Js), nbins)
for (j, J) in enumerate(Js)
    print(j, "-")
    # initialize model 
    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)
    H = DiagonalizedHamiltonian(m)
    E = H.E # energy eigenvalues 

    # calculate density of states 
    dos[j, :] = hist_counts(E, nbins=nbins)

    # calculate the IPR 
    ipr = IPR_real(H)
    # bin the results 
    ipr_bin = bin_results(ipr, nbins=nbins)
    iprs_real[j, :] = ipr_bin

    # calculate the IPR 
    ipr = IPR_momentum(H)
    # bin the results 
    ipr_bin = bin_results(ipr, nbins=nbins)
    iprs_k[j, :] = ipr_bin
end

nx = length(Js)
ny = nbins
hval = copy(iprs)

function get_colour(val; max_val)
    c = floor(Int, val / max_val * 100)
    if c == 0
        c = 1
    end
    cmap = cgrad(:magma, 100, categorical=true)
    return cmap[c]
end

p = plot(xlims=(minimum(Js), maximum(Js)), ylims=(0, ny + 1), grid=false, xaxis=:log10)
for (x, J) in enumerate(Js) # potential strength  
    for y in 1:ny # eigenstates 
        dos_xy = dos[x, y]
        val_xy = hval[x, y]
        # phase 
        c = get_colour(val_xy, max_val=maximum(hval))
        # onsite term  
        scatter!(p, [J], [y], ms=0.1 * abs(dos_xy), c=c, legend=:false, xaxis=:log10)
    end
end

# title!("⟨r⟩ for $L×$L×$L lattice")
plot(Js, iprs_real[:, 20], color="red", label=nothing)
plot!(Js, iprs_k[:, 20], color="blue", label=nothing)
scatter!(Js, iprs_real[:, 20], color="red", label="real space")
scatter!(Js, iprs_k[:, 20], color="blue", label="momentum space")


title!("IPR for $L×$L×$L lattice")