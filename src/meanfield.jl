using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations
using LoopVectorization
using ITensors

include("../src/results.jl")
include("../src/model.jl")

function pairfield_correlation(T; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, ϕ::Real=0, μ::Real, V0::Real, V1::Real=0, periodic::Bool=true, symmetry::String="s-wave")
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕ=ϕ, periodic=periodic)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    if symmetry == "s-wave" || symmetry == "s"
        M = swave(T, E=E, U=U, V0=V0)
    elseif symmetry == "d-wave" || symmetry == "d"
        M = dwave(T, L=L, E=E, U=U, V0=V0, V1=V1)
    else
        @error "Pairing symmetry not recognized"
        return
    end

    λ, Δ = calculate_λ_Δ(M)

    return λ, Δ
end

function decomposition_M(M)
    decomp, _ = partialschur(M, nev=1, tol=1e-6, which=LR())
    return decomp.R
end

function swave(T::Real; E, U, V0)
    println("s-wave configuration")
    N = size(U)[1]
    χ = zeros(N, N)

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
    PUU = reshape(PUU, N, N * N)
    UU = reshape(UU, N * N, N)
    χ = PUU * UU

    return -V0 * χ
end

function dwave_blocks(b_sites, d_sites; P, U, Uconj, V::Real, N::Int)
    Uconjb = [Uconj[r, :] for r in b_sites]
    Ud = [U[r, :] for r in d_sites]

    Uconjb = transpose(hcat(Uconjb...))
    Ud = transpose(hcat(Ud...))

    @tullio PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconjb[r, m]
    PUU = reshape(PUU, N, N * N)

    @tullio UU1[m, n, rprime] := U[rprime, n] * Ud[rprime, m]
    UU1 = reshape(UU1, N * N, N)
    χ1 = PUU * UU1

    @tullio UU2[m, n, rprime] := Ud[rprime, n] * U[rprime, m]
    UU2 = reshape(UU2, N * N, N)
    χ2 = PUU * UU2

    return -V / 2 .* (χ1 .+ χ2)
end

function dwave(T::Real; L, E, U, V0, V1)
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
        @show (b, d)

        # don't do the s-wave component (we've done it already)
        if !(b == 5 && d == 5)
            b_sites, d_sites = sites[b], sites[d]

            if b == 5 # only δ=0 term gets V0, not δ'=0! 
                Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V0, N=N)
            else # bond terms have potential V1 
                Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V1, N=N)
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
    # perform the decomposition 
    decomp, _ = partialschur(Hermitian(M), nev=1, tol=1e-6, which=LR())

    # extract the maximum eigenvector/value pair 
    maxev = decomp.Q[:, 1]
    λ = decomp.R[1]

    @show λ

    return λ, maxev
end

# function dwave_hermitian(T::Real; L, E, U, V0, V1)
#     println("d-wave configuration")
#     N = L^2
#     Uconj = conj.(U) # This is U*

#     # Initialize the M matrix 
#     M = Matrix{Matrix{Float64}}(undef, 5, 5)

#     # make the prefactor (1-fₙ-fₘ)/(Eₙ + Eₘ)
#     fs = fermi.(E, T)
#     fnm = zeros(N, N)
#     Enm = zeros(N, N)
#     for i in 1:N
#         fnm[:, i] = fs .+ fs[i]
#         Enm[:, i] = E .+ E[i]
#     end
#     P = (1 .- fnm) ./ Enm

#     # the s-wave sector. This is in the (5,5) block of M matrix
#     @tullio PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]
#     @tullio UU[m, n, rprime] := U[rprime, n] * U[rprime, m]
#     PUU = reshape(PUU, N, N * N)
#     UU = reshape(UU, N * N, N)
#     M[5, 5] = -V0 * PUU * UU

#     # make lists of the nearest-neighbour sites 
#     Rsites, Usites, Lsites, Dsites, onsites = [], [], [], [], []
#     for r in 1:N
#         nnr = [nearest_neighbours(r, L=L)...] # get the nearest neighbours
#         push!(Rsites, nnr[1])
#         push!(Usites, nnr[2])
#         push!(Lsites, nnr[3])
#         push!(Dsites, nnr[4])
#         push!(onsites, r)
#     end
#     sites = [Rsites, Usites, Lsites, Dsites, onsites]

#     # iterate through each of the 5×5 blocks
#     for bd in CartesianIndices(M)
#         (b, d) = Tuple(bd)
#         @show (b, d)

#         # bc matrix is Hermitian, we only have to fill in lower diagonal
#         # also, don't do the s-wave component (we've done it already)
#         if d <= b && !(b == 5 && d == 5)
#             b_sites, d_sites = sites[b], sites[d]

#             if b == 5 || d == 5 # the on-site terms get potential V=V0
#                 @error "IT MIGHT JUST BE b==5, not both"
#                 Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V0, N=N)
#             else # bond terms have potential V1 
#                 Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V1, N=N)
#             end

#             # fill in the matrix 
#             M[b, d] = Mblock
#             # if off-diagonal, fill in the hermitian conjugate block
#             if d != b
#                 M[d, b] = Mblock'
#             end
#         end
#     end
#     M = mortar(M)

#     return M
# end

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