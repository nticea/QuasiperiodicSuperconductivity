## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra, Arpack, Plots
using Profile
using ProgressBars
using CSV
using DataFrames

include("../../src/stiffness.jl")
include("../../src/meanfield.jl")


## PARAMETERS ##
L = 7 # the full system is L × L 
t = 1 # hopping 
Q = (√5 - 1) / 2
μ = 0
θ = π / 7
ϕx = 0.1234
ϕy = 0.5678
V0 = 1
V1 = -1.5
periodic = true
niter = 500
tol = 1e-15
T = 0

# get the BdG coefficients 
if V1 == 0
    U, V, E = BdG_coefficients_swave(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise, Δ_init=Δ_init)
else
    U, V, E = BdG_coefficients_dwave(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise, Δ_init=Δ_init)
end

# number of sites (L × L)
N, _ = size(U)

# the fermi distribution (vector over all sites)
f = fermi.(E; T=T)

# make lists of the nearest-neighbour sites 
Rsites, Usites, Lsites, Dsites, onsites = [], [], [], [], []
xcoords, ycoords = [], []
for r in 1:N
    nnr = [nearest_neighbours(r, L=L)...] # get the nearest neighbours
    push!(Rsites, nnr[1])
    push!(Usites, nnr[2])
    push!(Lsites, nnr[3])
    push!(Dsites, nnr[4])
    push!(onsites, r)

    x, y = site_to_coordinate(r, L=L)
    push!(xcoords, x)
    push!(ycoords, y)
end

sites = [Rsites, Usites, Lsites, Dsites, onsites]
coords = [xcoords, ycoords]

# the minimum q I can consider is 1/L
N, _ = size(U)
L = √N
qs = 2π / L * collect(1:npts)
Πs = zeros(npts, 4)

# the diagonals have ΔE=0, which causes divergence 
@einsimd En1n2[n1, n2] := E[n1] - E[n2] + δ * 1im
indices = findall(isequal(0), En1n2)
En1n2[indices] .= Inf

# get the fermi energies
@einsimd fn1n2[n1, n2] := f[n1] - f[n2]
@einsimd fn1n2_minus = copy(fn1n2) .- 1

# compute the prefactor
@einsimd pf_ad[n1, n2] := fn1n2[n1, n2] / En1n2[n1, n2]
@einsimd pf_bc[n1, n2] := fn1n2_minus[n1, n2] / En1n2[n1, n2]

# perform extrapolation q → 0
Threads.@threads for (i, q) in collect(enumerate(qs))
    Πs[i, :] = real.(Πq_debug(sites, coords; U=U, V=V, pf_ad=pf_ad, pf_bc=pf_bc, t=t, q=q))
end

function Πq_debug(sites, coords; U, V, pf_ab, pf_bc, t, q)
    N, _ = size(U)
    i_sites = sites[5] # these are the on-sites 

    Π = zeros(4)
    Threads.@threads for j in 1:2#1:4
        j_sites = sites[j]

        if j == 1 || j == 3
            # For Πxx, we set qx = 0 and take qy → 0
            qx = 0
            qy = q
        else
            # For Πyy, we set qy = 0 and take qx → 0
            qx = q
            qy = 0
        end

        a = Aq_debug((qx, qy), i_sites, j_sites, coords, U=U, V=V)
        b = Bq_debug((qx, qy), i_sites, j_sites, coords, U=U, V=V)
        c = Cq_debug((qx, qy), i_sites, j_sites, coords, U=U, V=V)
        d = Dq_debug((qx, qy), i_sites, j_sites, coords, U=U, V=V)

        Πxx = (a[n, m] + d[n, m]) * pf_ad[n, m] + (b[n, m] + c[n, m]) * pf_bc[n, m]

        Πxx = real.(Πxx)
        Π[j] = Πxx
    end
    Π[3] = copy(Π[1])
    Π[4] = copy(Π[2])
end

function Aq_debug(q, i_sites, j_sites, coords; U, V)
    # unpack some of the arguments 
    qx, qy = q
    x, y = coords

    # make the blocks of i sites and j sites 
    Uconj_i, U_i, Uconj_j, U_j = ij_blocks(i_sites, j_sites, U, conj.(U))

    # make the exponential prefactor 
    exp_pf = exp.(-1im .* (qx .* x + qy .* y))

    @einsimd A[n1, n2] := exp_pf[r] * (Uconj_j[r, n1] * U_i[r, n2] - Uconj_i[r, n1] * U_j[r, n2])

    return A
end