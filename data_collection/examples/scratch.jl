## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

# Parameters 
L = 7
J = 100
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
ϕx = 0
ϕy = 0
ϕz = 0
periodic = false
ndims = 3

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims)

function plot_eigenstates(m::ModelParams; slice::Int=1, E₀::Real=0)
    E, U = diagonalize_hamiltonian(m)
    potmat = real.(U[:, floor(Int, length(E) / 2)])

    if ndims == 2
        potmat = reshape(potmat, L, L)
    elseif ndims == 3
        potmat = reshape(potmat, L, L, L)
        potmat = potmat[:, :, slice]
    end

    numpts = 10
    cm = cgrad(:bwr, 2 * numpts + 1, categorical=true)

    function colour_gradient(x1::Int, x2::Int; arr)
        val = arr[x1, x2]
        max = 1
        idx = floor(Int, val / max * numpts + numpts + 1)
        return cm[idx]
    end

    plot(xaxis=" ", yaxis=" ", legend=false)
    heatmap([-1 1; 1 1], color=cm, visible=false)  # Dummy heatmap to generate colorbar

    # Create the colorbar

    h = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
    for x in 1:L
        for y in 1:L
            # onsite dot 
            scatter!(h, [x], [y], ms=10, c=colour_gradient(x, y, arr=potmat), legend=:false, aspect_ratio=:equal)
        end
    end

    xticks!(h, collect(1:2:L))
    yticks!(h, collect(1:2:L))
    xlabel!(h, "Site (x)")
    ylabel!(h, "Site (y)")

    return h
end

hs = []
for s in 1:L
    push!(hs, plot_eigenstates(m, slice=s, E₀=0))
end