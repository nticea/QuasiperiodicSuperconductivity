using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations
using LoopVectorization
using Polynomials
using FFTW

include("results.jl")
include("model.jl")
include("results.jl")

function pairfield_correlation(m::ModelParams; T::Real,
    checkpointpath::Union{String,Nothing}=nothing)

    E, U = diagonalize_hamiltonian(m, loadpath=checkpointpath)

    # s-wave case is faster 
    if V1 == 0
        M = swave(m, T, E=E, U=U)
        return calculate_λ_Δ(M)
    end

    M = dwave(m, T, E=E, U=U)
    λ, Δ = calculate_λ_Δ(M)

    # We've only explicitly calculated half the bonds  
    if ndims == 2
        Δ̃ = to_5N_LGE_Δ(Δ; L=m.L)
    elseif ndims == 3
        Δ̃ = to_7N_LGE_Δ(Δ; L=m.L)
    else
        println("$ndims dimensions not yet implemented")
    end

    return λ, Δ̃
end

function pairfield_susceptibility(m;
    T, symmetry::String, Λ::Union{Nothing,Real}=nothing,
    checkpointpath::Union{String,Nothing}=nothing)

    @warn "I need to benchmark this again after implementing 3D system sizes"

    E, U = diagonalize_hamiltonian(m, loadpath=checkpointpath)
    ndims = m.ndims

    if !isnothing(Λ)
        @warn "Keeping only states close to εf"
        sortidx = sortperm(E)
        Ẽ = E[sortidx]
        Ũ = U[:, sortidx]
        E = Ẽ[Ẽ.<Λ.&&Ẽ.>-Λ]
        U = Ũ[:, Ẽ.<Λ.&&Ẽ.>-Λ]
    end

    # Construct the pairfield susceptibility
    if symmetry == "s-wave"
        χ = swave_χ(m, T, E=E, U=U)
        χ0 = sum(χ)
    elseif symmetry == "d-wave"
        χ = dwave_χ(m, T, E=E, U=U)
        if ndims == 2
            N = L * L
            nblocks = 3
        elseif ndims == 3
            N = L * L * L
            nblocks = 4
        else
            println("$dims dimensions not implemented")
        end

        χ0 = zeros(nblocks, nblocks)
        for (b1, i) in enumerate(1:N:nblocks*N)
            for (b2, j) in enumerate(1:N:nblocks*N)
                χδδ = χ[i:(i+N-1), j:(j+N-1)]
                # take the sum of this guy 
                χ0[b1, b2] = sum(χδδ)
            end
        end
    else
        @error "Symmetry $symmetry not recognized"
    end

    return χ0
end

function uniform_susceptibility(m;
    T, symmetry::String="d-wave", Λ::Union{Nothing,Real}=nothing,
    checkpointpath::Union{String,Nothing}=nothing)

    ndims = m.ndims
    if ndims == 2
        N = L * L
        nblocks = 3
        if symmetry == "s-wave"
            δδp = [(1, 1)]
        else
            δδp = [(1, 1), (2, 2), (3, 3), (2, 3), (3, 2)]
        end
    elseif ndims == 3
        N = L * L * L
        nblocks = 4
        if symmetry == "s-wave"
            δδp = [(1, 1)]
        else
            δδp = [(1, 1), (2, 2), (3, 3), (4, 4), (2, 3), (3, 2),
                (2, 4), (4, 2), (3, 4), (4, 3)]
        end
    else
        println("$dims dimensions not implemented")
        return
    end

    E, U = diagonalize_hamiltonian(m, loadpath=checkpointpath)
    rvec = collect(1:N)

    # We need to transform each of the eigenvectors into 2D space! 
    function fourier_transform_U(u; minus=false)
        # first, map back to 2D space 
        if ndims == 2
            u = reshape(u, L, L)
        elseif ndims == 3
            u = reshape(u, L, L, L)
        else
            println("$dims dimensions not implemented")
            return
        end

        # take FT along all spatial dimensions 
        if !minus # this is the expression for U_{q}
            uq = fft(u)
        else # this is the expression for U_{-q}
            uq = conj.(fft(conj.(u)))
        end

        # reshape it back and normalize. FFTW does not normalize!!
        uq = reshape(uq, N) ./ √N

        return uq
    end

    # perform a Fourier transform along the real space dim
    Uq = fourier_transform_U.(eachcol(U))
    Uq = hcat(Uq...)
    Uminusq = fourier_transform_U.(eachcol(U), minus=true)
    Uminusq = hcat(Uminusq...)

    if !isnothing(Λ)
        @warn "Keeping only states close to εf"
        sortidx = sortperm(E)
        Ẽ = E[sortidx]
        Ũq = Uq[:, sortidx]
        Ũminusq = Uminusq[:, sortidx]
        E = Ẽ[Ẽ.<Λ.&&Ẽ.>-Λ]
        Uq = Ũq[:, Ẽ.<Λ.&&Ẽ.>-Λ]
        Uminusq = Ũminusq[:, Ẽ.<Λ.&&Ẽ.>-Λ]
    end

    Uq_conj = conj.(Uq)
    Uminusq_conj = conj.(Uminusq)

    # number of states we are keeping 
    N = size(Uq)[2]

    # make the prefactor
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    Pnm = (1 .- fnm) ./ Enm

    # I need the x and y components of q for the d-wave prefactor
    pfs = susceptibility_dwave_prefactors.(rvec, m=m)
    pfs = hcat(pfs...) # dimensions [δ] x [q]
    pfsneg = conj.(pfs)

    # multiply by prefactors 
    χ0 = zeros(nblocks, nblocks)
    for (δ, δp) in δδp
        # multiply with the prefactor 
        Uminusq_conj_δ = Uminusq_conj .* pfs[δ, :]
        Uminusq_δ = Uminusq .* pfsneg[δp, :]
        # create the terms 
        Tq = transpose(Uminusq_conj_δ) * Uq_conj
        Tl1 = transpose(Uminusq_δ) * Uq
        Tl2 = transpose(Uq) * Uminusq_δ
        # sum them together 
        Tl = Tl1 + Tl2
        @einsimd χ := 1 / (2 * L * L) * Pnm[n, m] * Tq[n, m] * Tl[n, m]
        # store the data 
        χ0[δ, δp] = real.(χ)
    end

    @show χ0
    return χ0
end

function susceptibility_dwave_prefactors(r::Int; m::ModelParams)
    qs = momentum_components(m, r=r)
    pfs = [exp(1im * q) for q in qs]
    pushfirst!(pfs, 1)
    return pfs
end

function decomposition_M(M)
    decomp, _ = partialschur(M, nev=1, tol=1e-6, which=LR())
    return decomp.R
end

function swave(m::ModelParams, T::Real; E, U)
    V0 = m.V0

    println("Computing s-wave configuration")
    Ntot, N = size(U)
    Uconj = conj.(U)

    # make the prefactor
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    P = (1 .- fnm) ./ Enm

    @einsimd PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]
    @einsimd UU[m, n, rprime] := U[rprime, n] * U[rprime, m]
    PUU = reshape(PUU, Ntot, N * N)
    UU = reshape(UU, N * N, Ntot)
    χ = PUU * UU

    return -V0 * χ
end

function swave_χ(m::ModelParams, T::Real; E, U)
    m_copy = copy(m)
    m_copy.V0 = -1
    return swave(m_copy, T, E=E, U=U)
end

function dwave_χ(m::ModelParams, T::Real; E, U)
    m_copy = copy(m)
    m_copy.V0 = -1
    m_copy.V1 = -1
    return dwave(m_copy, T, E=E, U=U)
end

function dwave_blocks(b_sites, d_sites; P, U, Uconj, V::Real, N::Int, Ntot::Int)
    Uconjb = [Uconj[r, :] for r in b_sites]
    Ud = [U[r, :] for r in d_sites]

    Uconjb = transpose(hcat(Uconjb...))
    Ud = transpose(hcat(Ud...))

    @tullio PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconjb[r, m]
    PUU = reshape(PUU, Ntot, N * N)

    @tullio UU1[m, n, rprime] := U[rprime, n] * Ud[rprime, m]
    UU1 = reshape(UU1, N * N, Ntot)
    χ1 = PUU * UU1

    @tullio UU2[m, n, rprime] := Ud[rprime, n] * U[rprime, m]
    UU2 = reshape(UU2, N * N, Ntot)
    χ2 = PUU * UU2

    return -V / 2 .* (χ1 .+ χ2)
end

function dwave(m::ModelParams, T::Real; E, U)
    V0, V1, ndims = m.V0, m.V1, m.ndims
    println("Computing d-wave configuration")
    Ntot, N = size(U)
    Uconj = conj.(U) # This is U^*

    # Initialize the M matrix 
    if ndims == 2
        M = Matrix{Matrix{Float64}}(undef, 3, 3)
    elseif ndims == 3
        M = Matrix{Matrix{Float64}}(undef, 4, 4)
    else
        println("$ndims dimensions not supported")
    end

    # make the prefactor (1-fₙ-fₘ)/(Eₙ + Eₘ)
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    P = (1 .- fnm) ./ Enm

    # the s-wave sector. This is in the (1,1) block of M matrix
    @tullio PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]
    @tullio UU[m, n, rprime] := U[rprime, n] * U[rprime, m]
    PUU = reshape(PUU, Ntot, N * N)
    UU = reshape(UU, N * N, Ntot)
    M[1, 1] = -V0 * PUU * UU

    # make lists of the nearest-neighbour sites 
    xsites, ysites, zsites, onsites = [], [], [], []
    for r in 1:Ntot
        nnr = [nearest_neighbours(r, m=m)...] # get the nearest neighbours
        push!(xsites, nnr[1])
        push!(ysites, nnr[2])
        push!(onsites, r)
        if ndims == 3
            push!(zsites, nnr[5])
        end
    end
    if ndims == 2
        sites = [onsites, xsites, ysites]
    elseif ndims == 3
        sites = [onsites, xsites, ysites, zsites]
    else
        println("$ndims dimensions is not supported")
    end

    # iterate through each of the blocks
    for bd in CartesianIndices(M)
        (b, d) = Tuple(bd)

        # don't do the s-wave component (we've done it already)
        if !(b == 1 && d == 1)
            b_sites, d_sites = sites[b], sites[d]

            if b == 1 # only δ=0 term gets V0, not δ'=0! 
                Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V0, N=N, Ntot=Ntot)
            else # bond terms have potential V1 
                # multiply by 2 bc we are only considering 1/2 of each direction
                Mblock = 2 * dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V1, N=N, Ntot=Ntot)
            end

            # fill in the matrix 
            M[b, d] = Mblock
        end
    end

    M = mortar(M)

    return M
end

function mortar(M::Matrix)
    (n, _) = size(M)
    (N, _) = size(M[1, 1])
    new_M = zeros(n * N, n * N)
    for b in 1:n
        for d in 1:n
            b_idx_start = (b - 1) * N + 1
            b_idx_end = b * N
            d_idx_start = (d - 1) * N + 1
            d_idx_end = d * N

            new_M[b_idx_start:b_idx_end, d_idx_start:d_idx_end] .= M[b, d]
        end
    end
    return new_M
end

function calculate_λ_Δ(M)
    try
        # perform the decomposition 
        decomp, _ = partialschur(M, nev=1, tol=1e-6, which=LR())

        # extract the maximum eigenvector/value pair 
        maxev = decomp.Q[:, 1]
        λ = decomp.R[1]

        @show λ
        return λ, maxev
    catch e
        @show e
        println("Trying a different decomposition")
        E, V = eigen(M)
        E = real.(E)
        sortidx = sortperm(E)
        E, V = E[sortidx], V[:, sortidx]

        λ = E[end]
        maxev = V[:, end]

        @show λ
        return λ, maxev
    end
end

function LGE_find_Tc(m; npts=5, tol=1e-4, L̃::Int=11)
    # find the min and max values based on the Tc of a smaller system
    m_small = ModelParams(L̃, m.t, m.Q, m.μ, m.θ, m.ϕx, m.ϕy, m.ϕz, m.V0, m.V1, m.J, m.periodic, m.ndims)
    Tc0, λ0, Δ0 = _LGE_find_Tc(m_small, npts=npts, tol=tol)
    if isnan(Tc0)
        println("No soln for small system size")
        return _LGE_find_Tc(m, min=0, max=1, npts=npts, tol=tol)
    end

    min = Tc0 - 0.2 * Tc0
    max = Tc0 + 0.4 * Tc0
    return _LGE_find_Tc(m, min=min, max=max, npts=npts, tol=tol)
end

function _LGE_find_Tc(m::ModelParams; min=0, max=1, npts=5, tol=1e-4, niter=10)
    λ0, Δ0 = pairfield_correlation(m, T=0)
    if λ0 < 1
        return NaN, λ0, Δ0
    end

    # else, we know that a value for Tc must exist 
    Tc = NaN
    for n in 1:niter
        Ts = LinRange(min, max, npts)
        λs = []
        for T in Ts
            # find λ at this temperature 
            λ, Δ = pairfield_correlation(m, T=T)

            @show T, λ
            # If λ is close enough to 1, return T as Tc 
            if abs.(λ - 1) < tol
                return T, λ, Δ
            end

            # if not, keep iterating 
            push!(λs, λ)

            if λ < 1
                break
            end
        end

        # find Tc 
        idxlist = sortperm(λs)
        knots = λs[idxlist]

        Tc = NaN
        @show length(knots)
        if length(knots) > 1
            Interpolations.deduplicate_knots!(knots, move_knots=true)
            try
                interp_linear = linear_interpolation(knots, Ts[idxlist])
                Tc = interp_linear(1)
            catch e
                if minimum(λs) > 1
                    max = 2 * max
                elseif maximum(λs) < 1
                    min = min / 2
                else
                    return NaN, λ0, Δ0
                end

                continue
            end
        else
            if minimum(λs) > 1
                max = 2 * max
            elseif maximum(λs) < 1
                min = min / 2
            else
                return NaN, λ0, Δ0
            end

            continue
        end

        # find λ at this temperature 
        λ, Δ = pairfield_correlation(m, T=Tc)

        # If λ is close enough to 1, return T as Tc 
        if abs.(λ - 1) < tol
            return Tc, λ, Δ
        end

        # else, keep iterating 
        min = Tc - 0.3 / n * Tc
        max = Tc + 0.3 / n * Tc
    end

    return Tc, λ0, Δ0
end