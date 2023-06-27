include("BdG.jl")
include("BdG_dwave.jl")

using Interpolations

function superfluid_stiffness(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0)

    # get the BdG coefficients 
    if V1 == 0
        U, V, E = BdG_coefficients_swave(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise)
    else
        U, V, E = BdG_coefficients_dwave(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise)
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

    K = kinetic_term(sites, U=U, V=V, E=E, f=f, t=t)
    Π = current_current_term(sites, coords, U=U, V=V, E=E, f=f)

    return K, Π

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

    @error "I might be missing a factor of 2 from the spin σ"
    @error "Check that I should be taking only positive eigenvalues"

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

        Kx = -t / N * (t1 + t2 + t3 + t4)

        push!(K, Kx)
    end

    return K
end

function current_current_term(sites, coords; U, V, E, f)
    N, _ = size(U)
    qtol = 0
    @error "Need to do a finite size extrapolation"

    i_sites = sites[5] # these are the on-sites 

    # the diagonals have ΔE=0, which causes divergence 
    @einsimd En1n2[n1, n2] := (E[n1] - E[n2])
    En1n2[diagind(En1n2)] .= Inf

    Π = []
    for j_sites in sites[1:4]
        A = Aq((qtol, qtol), i_sites, j_sites, coords, U=U)
        D = Dq((-qtol, -qtol), i_sites, j_sites, coords, V=V)

        Aconj = conj.(A)

        @einsum Πxx := 1 / N * (A[n1, n2] * (Aconj[n1, n2] + D[n1, n2]) * (f[n1] - f[n2]) / En1n2[n1, n2])

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

    # make the term 
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

    # make the term 
    exp_pf = exp.(-1im .* (qx .* x + qy .* y))
    @einsimd D[n1, n2] := exp_pf[r] * (V_j[r, n1] * Vconj_i[r, n2] - V_i[r, n1] * Vconj_j[r, n2])

    return D
end

function fermi(ε; T)
    1 / (1 + exp(ε / T))
end