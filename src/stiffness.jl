include("BdG.jl")
include("BdG_dwave.jl")

using Interpolations
using Polynomials

function superfluid_stiffness_finiteT(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0, Δ_init)

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

    K = @time kinetic_term(sites, U=U, V=V, E=E, f=f, t=t)
    Π = current_current_term(sites, coords, U=U, V=V, E=E, f=f, t=t)

    return K, Π
end

function superfluid_stiffness(; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0, npts::Int=5, Δ_init)

    Ts = expspace(-1, -9, npts)
    Ds = zeros(npts, 4)
    # collect data points at various T 
    for (i, T) in enumerate(Ts)
        print(i, "-")
        K, Π = superfluid_stiffness_finiteT(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise, Δ_init=Δ_init)
        @show K, Π
        @show -K + Π
        Ds[i, :] = -K + Π
    end

    # perform extrapolation T → 0
    Ds_extrapolated = []
    for x in 1:4
        model = Polynomials.fit(Ts, Ds[:, x], npts)
        push!(Ds_extrapolated, model(0))
    end

    return Ds_extrapolated
end

function ij_blocks(i_sites, j_sites, A, Aconj)
    Aconj_i = [Aconj[r, :] for r in i_sites]
    A_i = [A[r, :] for r in i_sites]
    Aconj_j = [Aconj[r, :] for r in j_sites]
    A_j = [A[r, :] for r in j_sites]

    Aconj_i = transpose(hcat(Aconj_i...))
    A_i = transpose(hcat(A_i...))
    Aconj_j = transpose(hcat(Aconj_j...))
    A_j = transpose(hcat(A_j...))

    return Aconj_i, A_i, Aconj_j, A_j
end

function kinetic_term(sites; U, V, E, f, t)

    N, _ = size(U)
    i_sites = sites[5] # these are the on-sites 

    # take only the positive eigenvalues
    En_idx = findall(x -> x >= 0, E)
    U = U[:, En_idx]
    V = V[:, En_idx]
    f = f[En_idx]

    K = []
    for j_sites in sites[1:4] #iterate through the nearest neighbours in all directions 
        # make the blocks of i sites and j sites 
        Uconj_i, U_i, Uconj_j, U_j = ij_blocks(i_sites, j_sites, U, conj.(U))
        Vconj_i, V_i, Vconj_j, V_j = ij_blocks(i_sites, j_sites, V, conj.(V))

        # calculate the terms 
        @einsimd t1 := f[n] * Uconj_j[r, n] * U_i[r, n]
        @einsimd t2 := f[n] * Uconj_i[r, n] * U_j[r, n]
        @einsimd t3 := (1 .- f[n]) * V_j[r, n] * Vconj_i[r, n]
        @einsimd t4 := (1 .- f[n]) * V_i[r, n] * Vconj_j[r, n]

        Kx = -2 * t / N * (t1 + t2 + t3 + t4) # factor of 2 from spin 

        push!(K, Kx)
    end

    return K
end

function current_current_term(sites, coords; U, V, E, f, t, npts=5, δ=1e-8)
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

    # compute the prefactor
    @einsimd pf[n1, n2] := fn1n2[n1, n2] / En1n2[n1, n2]

    # perform extrapolation q → 0
    for (i, q) in enumerate(qs)
        print(i, "-")
        Πs[i, :] = @time real.(Πq(sites, coords; U=U, V=V, pf=pf, t=t, q=q))
    end

    Πs_extrapolated = []
    for x in 1:4
        Πxx = Πs[:, x]
        model = Polynomials.fit(qs, Πxx, npts)
        push!(Πs_extrapolated, model(0))
    end

    return Πs_extrapolated
end

function Πq(sites, coords; U, V, pf, t, q)
    N, _ = size(U)
    i_sites = sites[5] # these are the on-sites 

    Π = []
    for j in 1:4
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

        A = Aq((qx, qy), i_sites, j_sites, coords, U=U)
        D = Dq((-qx, -qy), i_sites, j_sites, coords, V=V)
        Aconj = conj.(A)

        # @einsimd Πxx := 2 * t^2 / N * (A[n1, n2] * (Aconj[n1, n2] + D[n1, n2]) * (fn1n2[n1, n2]) / En1n2[n1, n2])
        AplusD = Aconj + D
        @einsimd Πxx := 2 * t^2 / N * A[n1, n2] * AplusD[n1, n2] * pf[n1, n2]
        Πxx = real.(Πxx)
        push!(Π, Πxx)
    end

    return Π
end

function Aq(q, i_sites, j_sites, coords; U)
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

function Dq(q, i_sites, j_sites, coords; V)

    # unpack some of the arguments 
    qx, qy = q
    x, y = coords

    # make the blocks of i sites and j sites 
    Vconj_i, V_i, Vconj_j, V_j = ij_blocks(i_sites, j_sites, V, conj.(V))

    # make the exponential prefactor
    exp_pf = exp.(-1im .* (qx .* x + qy .* y))

    # factor of 2 is for the spin ??
    @einsimd D[n1, n2] := exp_pf[r] * (V_j[r, n1] * Vconj_i[r, n2] - V_i[r, n1] * Vconj_j[r, n2])

    return D
end

function fermi(ε; T)
    1 / (1 + exp(ε / T))
end