
## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots, Statistics, ExponentialUtilities
include("../../src/model.jl")

## PARAMETERS ## 
L = 15
Js = [0, 0.1, 0.5, 1, 5, 10]
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = -1
V1 = -1
ϕx = 0
ϕy = 0
ϕz = 0
periodic = false
disorder = false
ndims = 3
numsteps = 5000
ΔT = 1e-2

# construct the Hamiltonian 
m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=0, V1=0, J=0, periodic=periodic, ndims=ndims, disorder=disorder)
# make the initial state 
Lsite = L
if m.ndims == 2
    rLLL = coordinate_to_site(Lsite, Lsite, m=m)
elseif m.ndims == 3
    rLLL = coordinate_to_site(Lsite, Lsite, Lsite, m=m)
end
r₀ = zeros(numsites(m))
r₀[rLLL] = 1

# calculate all positions and differences 
function mindist(x, x̃, x0)
    dx1 = (x - x0)^2
    dx̃1 = (x̃ - x0)^2
    if periodic
        return minimum([dx1, dx̃1])
    else
        return dx1
    end
end

function distance(r; r0)
    x, y, z = site_to_coordinate(r, m=m)
    x0, y0, z0 = site_to_coordinate(r0, m=m)
    x̃, ỹ, z̃ = x + L, y + L, z + L
    dx2 = mindist(x, x̃, x0)
    dy2 = mindist(y, ỹ, y0)
    dz2 = mindist(z, z̃, z0)
    return dx2 + dy2 + dz2
end

dists = []
for r in collect(1:numsites(m))
    push!(dists, distance(r, r0=rLLL))
end

# evolve the state forward 
δrs = zeros(numsteps, length(Js))
for (Jᵢ, J) in enumerate(Js)
    # construct the Hamiltonian 
    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=0, V1=0, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
    H0 = noninteracting_hamiltonian(m, scale_model=false)

    # initialize the previous state as the initial state 
    r_prev = copy(r₀)
    for Tᵢ in 1:numsteps
        print("$Tᵢ-")

        rᵢ = expv(-1im * ΔT, H0, r_prev)
        r_prev = copy(rᵢ)

        # compute difference (these methods give same result)
        δrs[Tᵢ, Jᵢ] = sum(dists .* abs2.(rᵢ))
    end
end

if periodic
    BC = "PBC"
else
    BC = "OBC"
end
if disorder
    pot = "disorder"
else
    pot = "QP"
end

p = plot(xlabel="t", ylabel="δr²(t)", title="Transport (L=$L, $BC, $pot potential)", legend=:left)
Ts = collect(ΔT:ΔT:ΔT*numsteps)
for (Jᵢ, J) in enumerate(Js)
    p = plot!(p, Ts, δrs[:, Jᵢ], label="J=$J", xaxis=:log10, yaxis=:log10)
end


