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

include("../src/results.jl")
include("../src/model.jl")

function pairfield_correlation(T; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0, μ::Real, V0::Real, V1::Real=0, periodic::Bool=true)
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    if V1 == 0
        M = swave(T, E=E, U=U, V0=V0)
    else
        M = dwave(T, L=L, E=E, U=U, V0=V0, V1=V1)
    end

    λ, Δ = calculate_λ_Δ(M)

    # the result is only 3N. Need to make it 5N 
    Δ̃ = to_5N_LGE_Δ(Δ; L=L)

    return λ, Δ̃
end

function pairfield_susceptibility(T, symmetry::String; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0, μ::Real, periodic::Bool=true, Λ::Union{Nothing,Float64}=nothing)
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

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
        return swave_χ(T, E=E, U=U)
    elseif symmetry == "d-wave"
        return dwave_χ(T, L=L, E=E, U=U)
    else
        @error "Symmetry $symmetry not recognized"
    end
end

function uniform_susceptibility(T, symmetry::String; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0, μ::Real, periodic::Bool=true, Λ::Union{Nothing,Float64}=nothing)
    # Construct the non-interacting Hamiltonian matrix
    χ = pairfield_susceptibility(T, symmetry, L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic, Λ=Λ)
    if symmetry == "s-wave"
        @einsimd Δ := χ[r, r]
    end
end

function return_M(T; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0, μ::Real, V0::Real, V1::Real=0, periodic::Bool=true)
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    if V1 == 0
        M = swave(T, E=E, U=U, V0=V0)
    else
        M = dwave(T, L=L, E=E, U=U, V0=V0, V1=V1)
    end

    return M
end

function decomposition_M(M)
    decomp, _ = partialschur(M, nev=1, tol=1e-6, which=LR())
    return decomp.R
end

function swave_χ(T::Real; E, U)
    return swave(T, E=E, U=U, V0=-1)
end

function swave(T::Real; E, U, V0)
    println("s-wave configuration")
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

function dwave_χ(T::Real; L, E, U)
    return dwave(T, L=L, E=E, U=U, V0=-1, V1=-1)
end

function dwave(T::Real; L, E, U, V0, V1)
    println("d-wave configuration")
    N = size(U)[2]
    Ntot = L^2
    Uconj = conj.(U) # This is U^*

    # Initialize the M matrix 
    M = Matrix{Matrix{Float64}}(undef, 3, 3)

    # make the prefactor (1-fₙ-fₘ)/(Eₙ + Eₘ)
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    P = (1 .- fnm) ./ Enm

    # the s-wave sector. This is in the (5,5) block of M matrix
    @tullio PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]
    @tullio UU[m, n, rprime] := U[rprime, n] * U[rprime, m]
    PUU = reshape(PUU, Ntot, N * N)
    UU = reshape(UU, N * N, Ntot)
    M[3, 3] = -V0 * PUU * UU

    # make lists of the nearest-neighbour sites 
    # Rsites, Usites, Lsites, Dsites, onsites = [], [], [], [], []
    Rsites, Usites, onsites = [], [], []
    for r in 1:Ntot
        nnr = [nearest_neighbours(r, L=L)...] # get the nearest neighbours
        push!(Rsites, nnr[1])
        push!(Usites, nnr[2])
        # push!(Lsites, nnr[3])
        # push!(Dsites, nnr[4])
        push!(onsites, r)
    end
    # sites = [Rsites, Usites, Lsites, Dsites, onsites]
    sites = [Rsites, Usites, onsites]

    # iterate through each of the 3×3 blocks
    for bd in CartesianIndices(M)
        (b, d) = Tuple(bd)

        # don't do the s-wave component (we've done it already)
        if !(b == 3 && d == 3)
            b_sites, d_sites = sites[b], sites[d]

            if b == 3 # only δ=0 term gets V0, not δ'=0! 
                Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V0, N=N, Ntot=Ntot)
            else # bond terms have potential V1 
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

function dwave_hermitian(T::Real; L, E, U, V0, V1)
    println("d-wave configuration")
    N = L^2
    Uconj = conj.(U) # This is U*

    # Initialize the M matrix 
    M = Matrix{Matrix{Float64}}(undef, 5, 5)

    # make the prefactor (1-fₙ-fₘ)/(Eₙ + Eₘ)
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    P = (1 .- fnm) ./ Enm

    # the s-wave sector. This is in the (5,5) block of M matrix
    @tullio PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]
    @tullio UU[m, n, rprime] := U[rprime, n] * U[rprime, m]
    PUU = reshape(PUU, N, N * N)
    UU = reshape(UU, N * N, N)
    M[5, 5] = -V0 * PUU * UU

    # make lists of the nearest-neighbour sites 
    Rsites, Usites, Lsites, Dsites, onsites = [], [], [], [], []
    for r in 1:N
        nnr = [nearest_neighbours(r, L=L)...] # get the nearest neighbours
        push!(Rsites, nnr[1])
        push!(Usites, nnr[2])
        push!(Lsites, nnr[3])
        push!(Dsites, nnr[4])
        push!(onsites, r)
    end
    sites = [Rsites, Usites, Lsites, Dsites, onsites]

    # iterate through each of the 5×5 blocks
    for bd in CartesianIndices(M)
        (b, d) = Tuple(bd)

        # bc matrix is Hermitian, we only have to fill in lower diagonal
        # also, don't do the s-wave component (we've done it already)
        if d <= b && !(b == 5 && d == 5)
            b_sites, d_sites = sites[b], sites[d]

            if b == 5 || d == 5 # the on-site terms get potential V=V0
                @error "IT MIGHT JUST BE b==5, not both"
                Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V0, N=N)
            else # bond terms have potential V1 
                Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V1, N=N)
            end

            # fill in the matrix 
            M[b, d] = Mblock
            # if off-diagonal, fill in the hermitian conjugate block
            if d != b
                M[d, b] = Mblock'
            end
        end
    end
    M = mortar(M)

    return M
end

function LGE_spectrum(T; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0, μ::Real, V0::Real, V1::Real=0, periodic::Bool=true)
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    if V1 == 0
        M = swave(T, E=E, U=U, V0=V0)
    else
        M = dwave(T, L=L, E=E, U=U, V0=V0, V1=V1)
    end

    # diagonalize M fully
    vals, vecs = eigen(M)

    @assert maximum(imag.(vals)) < 1e-14
    vals = real.(vals)

    # sort the eigenvalues
    sort!(vals)

    return vals
end

function LGE_find_Tc(; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0, μ::Real, V0::Real, V1::Real=0, periodic::Bool=true, npts=5, tol=1e-4, niter=10, L̃::Int=11)
    # find the min and max values based on the Tc of a smaller system 
    Tc0, λ0, Δ0 = _LGE_find_Tc(L=L̃, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, npts=npts, tol=tol)
    if isnan(Tc0)
        println("No soln for small system size")
        return _LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, min=0, max=1, npts=npts, tol=tol)
    end

    min = Tc0 - 0.2 * Tc0
    max = Tc0 + 0.4 * Tc0
    return _LGE_find_Tc(L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic, min=min, max=max, npts=npts, tol=tol)
end

function _LGE_find_Tc(; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}, ϕx::Real, ϕy::Real, μ::Real, V0::Real, V1::Real, periodic::Bool=true, min=0, max=1, npts=5, tol=1e-4, niter=10)
    λ0, Δ0 = pairfield_correlation(0, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)
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
            λ, Δ = pairfield_correlation(T, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)

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
        λ, Δ = pairfield_correlation(Tc, L=L, t=t, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy, μ=μ, V0=V0, V1=V1, periodic=periodic)

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

# function dwave_old(T::Real; L, E, U, V0, V1)
#     println("d-wave configuration")
#     N = L^2

#     # make the prefactor
#     fs = fermi.(E, T)
#     fnm = zeros(N, N)
#     Enm = zeros(N, N)
#     for i in 1:N
#         fnm[:, i] = fs .+ fs[i]
#         Enm[:, i] = E .+ E[i]
#     end
#     P = (1 .- fnm) ./ Enm

#     Uconj = conj.(U)

#     # put together the M matrix 
#     M = zeros(N, 5, N, 5)
#     for r in 1:N
#         for rp in 1:N

#             # nearest neighbours of r 
#             nn = [nearest_neighbours(r, L=L)...]
#             push!(nn, r)

#             # nearest neighbours of r' 
#             nnrp = [nearest_neighbours(rp, L=L)...]
#             push!(nnrp, rp)

#             # iterate through δ
#             for (idx_δ, δ) in enumerate(nn)

#                 # iterate through δ'
#                 for (idx_δp, δp) in enumerate(nnrp)

#                     Ua = Uconj[r, :]
#                     Ub = Uconj[δ, :]
#                     Uc = U[rp, :]
#                     Ud = U[δp, :]

#                     @einsimd χ1 := P[n, m] * Ua[n] * Uc[n] * Ub[m] * Ud[m]
#                     @einsimd χ2 := P[n, m] * Ua[n] * Ud[n] * Ub[m] * Uc[m]

#                     χ = 1 / 2 * (χ1 + χ2)

#                     if idx_δ == 5 #|| idx_δp == 5
#                         M[r, idx_δ, rp, idx_δp] = -V0 * χ
#                     else
#                         M[r, idx_δ, rp, idx_δp] = -V1 * χ
#                     end
#                 end

#             end
#         end
#     end

#     i1 = Index(N)
#     i2 = Index(5)
#     i3 = Index(N)
#     i4 = Index(5)
#     M = ITensor(M, i1, i2, i3, i4)

#     C1 = combiner(i1, i2; tags="c1")
#     C2 = combiner(i3, i4; tags="c2")

#     M = M * C1 * C2

#     M = array(M)

#     # reshape M 
#     # M = reshape(M, 5 * N, 5 * N)

#     @show maximum(M - M')

#     return M
# end