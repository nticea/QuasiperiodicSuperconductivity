include("BdG.jl")

using Interpolations
using Polynomials

function superfluid_stiffness_finiteT(m::ModelParams; T::Real, niter::Int=100, tol::Real=1e-12, Δ_init=nothing)

    U, V, E, Δ = BdG_coefficients(m, T=T, niter=niter, tol=tol, Δ_init=Δ_init)

    # number of sites
    N, _ = size(U)
    ndims = m.ndims

    # the fermi distribution (vector over all sites)
    f = fermi.(E; T=T)

    # make lists of the nearest-neighbour sites 
    Rsites, Usites, zsites, onsites = [], [], [], []
    xcoords, ycoords, zcoords = [], [], []
    for r in 1:N
        nnr = [nearest_neighbours(r, m=m)...] # get the nearest neighbours
        push!(onsites, r)
        push!(Rsites, nnr[3])
        push!(Usites, nnr[2])
        push!(zsites, nnr[5])

        if ndims == 2
            x, y = site_to_coordinate(r, m=m)
            push!(xcoords, x)
            push!(ycoords, y)
        elseif ndims == 3
            x, y, z = site_to_coordinate(r, m=m)
            push!(xcoords, x)
            push!(ycoords, y)
            push!(zcoords, z)
        else
            println("sorry, bruv")
            return
        end
    end

    if ndims == 2
        sites = [onsites, Rsites, Usites]
        coords = [xcoords, ycoords]
    elseif ndims == 3
        sites = [onsites, Rsites, Usites, zsites]
        coords = [xcoords, ycoords, zcoords]
    else
        println("tough.")
        return
    end

    @time K = kinetic_term(sites, U=U, V=V, E=E, f=f, m=m)
    @time Π = current_current_term(sites, coords, U=U, V=V, E=E, f=f, m=m)

    return K, Π, Δ
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

function kinetic_term(sites; U, V, E, f, m)
    t, ndims = m.t, m.ndims

    N, _ = size(U)
    i_sites = sites[1] # these are the on-sites 

    # take only the positive eigenvalues
    En_idx = findall(x -> x >= 0, E)
    U = U[:, En_idx]
    V = V[:, En_idx]
    f = f[En_idx]

    if ndims == 2
        K = zeros(2)
    elseif ndims == 3
        K = zeros(3)
    else
        println("rough")
        return
    end

    for j in 2:length(sites) #iterate through the nearest neighbours in all directions 
        j_sites = sites[j]

        # make the blocks of i sites and j sites 
        Uconj_i, U_i, Uconj_j, U_j = ij_blocks(i_sites, j_sites, U, conj.(U))
        Vconj_i, V_i, Vconj_j, V_j = ij_blocks(i_sites, j_sites, V, conj.(V))

        fminus = 1 .- f
        # calculate the terms 
        @einsimd t1 := f[n] * Uconj_j[r, n] * U_i[r, n]
        @einsimd t2 := f[n] * Uconj_i[r, n] * U_j[r, n]
        @einsimd t3 := fminus[n] * V_j[r, n] * Vconj_i[r, n]
        @einsimd t4 := fminus[n] * V_i[r, n] * Vconj_j[r, n]

        Kx = -2 * t / N * (t1 + t2 + t3 + t4) # factor of 2 from spin 

        K[j-1] = Kx
    end

    return K
end

function current_current_term(sites, coords; U, V, E, f, m::ModelParams, npts=5, δ=1e-8, porder=3)
    ndims = m.ndims

    # the minimum q I can consider is 1/L
    qs = expspace(1, -9, npts)
    if ndims == 2
        Πs = zeros(npts, 2)
    elseif ndims == 3
        Πs = zeros(npts, 3)
    else
        println(":(")
        return
    end

    # the diagonals have ΔE=0, which causes divergence 
    @einsimd En1n2[n1, n2] := E[n1] - E[n2] + δ * 1im
    indices = findall(isequal(0), En1n2)
    En1n2[indices] .= Inf

    # get the fermi energies
    @einsimd fn1n2[n1, n2] := f[n1] - f[n2]

    # compute the prefactor
    @einsimd pf[n1, n2] := fn1n2[n1, n2] / En1n2[n1, n2]

    # perform extrapolation qy → 0
    Threads.@threads for (i, q) in collect(enumerate(qs))
        print(i, "-")
        Πs[i, :] = real.(Πq(sites, coords; U=U, V=V, pf=pf, q=q, m=m))
    end

    @show Πs[:, 1]

    if npts <= porder
        porder = npts - 2
    end
    Πs_extrapolated = []
    for x in 1:size(Πs)[2]
        Πxx = Πs[:, x]
        model = Polynomials.fit(qs, Πxx, porder)
        push!(Πs_extrapolated, model(0))
    end

    @show Πs_extrapolated[1]

    return Πs_extrapolated
end

function Πq(sites, coords; U, V, pf, q, m::ModelParams)
    t, ndims = m.t, m.ndims

    N, _ = size(U)
    i_sites = sites[1] # these are the on-sites 

    if ndims == 2
        Π = zeros(2)
    elseif ndims == 3
        Π = zeros(3)
    else
        println("sad.")
        return
    end

    Threads.@threads for j in 2:length(sites)
        j_sites = sites[j]

        if ndims == 2
            if j == 2
                # For Πxx, we set qx = 0 and take qy → 0
                qx = 0
                qy = q
            elseif j == 3
                # For Πyy, we set qy = 0 and take qx → 0
                qx = q
                qy = 0
            else
                println("Something is off with indexing...")
            end
            A = Aq((qx, qy), i_sites, j_sites, coords, U=U, m=m)
            D = Dq((-qx, -qy), i_sites, j_sites, coords, V=V, m=m)
        elseif ndims == 3
            if j == 2
                # For Πxx, we set qx = 0 and take qy → 0
                qx = 0
                qy = q
                qz = q
            elseif j == 3
                # For Πyy, we set qy = 0 and take qx → 0
                qx = q
                qy = 0
                qz = q
            elseif j == 4
                qx = q
                qy = q
                qz = 0
            else
                println("Not sure how you got here. Good job!")
                return
            end
            A = Aq((qx, qy, qz), i_sites, j_sites, coords, U=U, m=m)
            D = Dq((-qx, -qy, -qz), i_sites, j_sites, coords, V=V, m=m)
        else
            println("boop")
            return
        end
        
        Aconj = conj.(A)
        AplusD = Aconj + D
        @einsimd Πxx := 2 * t^2 / N * A[n1, n2] * AplusD[n1, n2] * pf[n1, n2]
        Πxx = real.(Πxx)
        Π[j-1] = Πxx
    end
    return Π
end

function Aq(q, i_sites, j_sites, coords; U, m::ModelParams)
    if m.ndims == 2
        # unpack some of the arguments 
        qx, qy = q
        x, y = coords
        # make the exponential prefactor 
        exp_pf = exp.(-1im .* (qx .* x + qy .* y))
    elseif m.ndims == 3
        # unpack some of the arguments 
        qx, qy, qz = q
        x, y, z = coords
        # make the exponential prefactor 
        exp_pf = exp.(-1im .* (qx .* x + qy .* y + qz .* z))
    else
        println("no.")
        return
    end

    # make the blocks of i sites and j sites 
    Uconj_i, U_i, Uconj_j, U_j = ij_blocks(i_sites, j_sites, U, conj.(U))

    # this is way faster! 
    U1 = Uconj_j .* exp_pf
    U2 = Uconj_i .* exp_pf
    t1 = transpose(U1) * U_i
    t2 = transpose(U2) * U_j
    A = t1 - t2

    return A
end

function Dq(q, i_sites, j_sites, coords; V, m::ModelParams)
    if m.ndims == 2
        # unpack some of the arguments 
        qx, qy = q
        x, y = coords
        # make the exponential prefactor 
        exp_pf = exp.(-1im .* (qx .* x + qy .* y))
    elseif m.ndims == 3
        # unpack some of the arguments 
        qx, qy, qz = q
        x, y, z = coords
        # make the exponential prefactor 
        exp_pf = exp.(-1im .* (qx .* x + qy .* y + qz .* z))
    else
        println("no.")
        return
    end
    # make the blocks of i sites and j sites 
    Vconj_i, V_i, Vconj_j, V_j = ij_blocks(i_sites, j_sites, V, conj.(V))

    # this is way faster! 
    U1 = V_j .* exp_pf
    U2 = V_i .* exp_pf
    t1 = transpose(U1) * Vconj_i
    t2 = transpose(U2) * Vconj_j
    D = t1 - t2

    return D
end

function fermi(ε; T)
    1 / (1 + exp(ε / T))
end