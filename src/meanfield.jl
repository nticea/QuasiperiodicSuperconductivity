using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations

include("../src/results.jl")
include("../src/model.jl")

function λmax(T; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, μ::Real, V0::Real, V1::Real=0, periodic::Bool=true, symmetry::String="s-wave")
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, periodic=periodic)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    if symmetry == "s-wave" || symmetry == "d"
        M = swave(T, E=E, U=U, V0=V0)
    elseif symmetry == "d-wave" || symmetry == "d"
        @time M = dwave(T, L=L, E=E, U=U, V0=V0, V1=V1)
    else
        @error "Pairing symmetry not recognized"
        return
    end

    # Calculate Tc by finding the eigenvalues of M
    @time λs = decomposition_M(M)

    ## TODO: DIAGONALIZE EACH BLOCK INDEPENDENTLY

    return λs[1]
end

function pairfield_correlation(T; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, μ::Real, V0::Real, V1::Real=0, periodic::Bool=true, symmetry::String="s-wave")
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, periodic=periodic)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    if symmetry == "s-wave" || symmetry == "d"
        M = swave(T, E=E, U=U, V0=V0)
    elseif symmetry == "d-wave" || symmetry == "d"
        @time M = dwave(T, L=L, E=E, U=U, V0=V0, V1=V1)
    else
        @error "Pairing symmetry not recognized"
        return
    end

    return M
end

function decomposition_M(M)
    decomp, _ = partialschur(M, nev=1, tol=1e-6, which=LM())
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

    return V0 * χ
end

## TODO:: CONSTRUCT EACH BLOCK INDEPENDENTLY
# The fast version of the code would like something like
# FOREACH BLOCK (that is to say, δ=[(0,0),(1,1),(1,0),(0,1),(-1,0),(0,-1)]):
# @einsimd PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]
# @einsimd UU[m, n, rprime] := U[rprime, n] * U[rprime+δ, m]
# PUU = reshape(PUU, N, N * N)
# UU = reshape(UU, N * N, N)

# function dwave(T::Real; L, E, U, V0, V1)
#     println("d-wave configuration")
#     N = L^2
#     Vs = [V1, V1, V1, V1, V0]
#     Uconj = conj.(U) # This is U*

#     # make the prefactor (1-fₙ-fₘ)/(Eₙ + Eₘ)
#     fs = fermi.(E, T)
#     fnm = zeros(N, N)
#     Enm = zeros(N, N)
#     for i in 1:N
#         fnm[:, i] = fs .+ fs[i]
#         Enm[:, i] = E .+ E[i]
#     end
#     P = (1 .- fnm) ./ Enm

#     # The M matrix has dimension 5L^2 × 5L^2 
#     M = zeros(N, 5, N, 5)

#     Threads.@threads for r in 1:N # iterate through sites r 

#         # get the nearest neighbours of r 
#         nnr = [nearest_neighbours(r, L=L)...]

#         # concatenate to this list the site itself (s-wave component)
#         push!(nnr, r)

#         # iterate through the nearest neighbours of r 
#         for (idx_rd, rd) in enumerate(nnr)

#             for rp in 1:N # iterate through sites r'

#                 # get the nearest neighbours
#                 nnrp = [nearest_neighbours(rp, L=L)...]

#                 # concatenate to this list the site itself (s-wave component)
#                 push!(nnrp, rp)

#                 # iterate through sites r'+δ'
#                 for (idx_rpd, rpd) in enumerate(nnrp)

#                     Ua = Uconj[r, :]
#                     Ub = Uconj[rd, :]
#                     Uc = U[rp, :]
#                     Ud = U[rpd, :]

#                     @einsimd Uterm1[n, m] := Ua[n] * Uc[n] * Ub[m] * Ud[m]
#                     @einsimd Uterm2[n, m] := Ua[n] * Ud[n] * Ub[m] * Uc[m]
#                     Uterm = Uterm1 .+ Uterm2

#                     @einsimd χ := P[n, m] * Uterm[n, m]

#                     if idx_rd == 5 || idx_rpd == 5
#                         V = V0
#                     else
#                         V = V1
#                     end

#                     M[r, idx_rd, rp, idx_rpd] = 1 / 2 * V * χ
#                 end
#             end
#         end
#     end

#     # reshape M 
#     M = reshape(M, 5 * N, 5 * N)

#     return M
# end

function dwave_blocks(b_sites, d_sites; P, U, Uconj, V::Real, N::Int)
    Uconjb = [Uconj[r, :] for r in b_sites]
    Ud = [U[r, :] for r in d_sites]

    Uconjb = transpose(hcat(Uconjb...))
    Ud = transpose(hcat(Ud...))

    @einsimd PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconjb[r, m]
    PUU = reshape(PUU, N, N * N)

    @einsimd UU[m, n, rprime] := U[rprime, n] * Ud[rprime, m]
    UU = reshape(UU, N * N, N)
    χ1 = PUU * UU

    @einsimd UU[m, n, rprime] := Ud[rprime, n] * U[rprime, m]
    UU = reshape(UU, N * N, N)
    χ2 = PUU * UU

    return V / 2 .* (χ1 .+ χ2)
end

function dwave(T::Real; L, E, U, V0, V1)
    println("d-wave configuration")
    N = L^2
    Uconj = conj.(U) # This is U*

    M = zeros(5 * N, 5 * N)

    # make the prefactor (1-fₙ-fₘ)/(Eₙ + Eₘ)
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    P = (1 .- fnm) ./ Enm

    # the s-wave sector 
    # @einsimd PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]
    # @einsimd UU[m, n, rprime] := U[rprime, n] * U[rprime, m]
    # PUU = reshape(PUU, N, N * N)
    # UU = reshape(UU, N * N, N)
    # M[4*N+1:end, 4*N+1:end] .= V0 * PUU * UU

    # the nearest-neighbour sites
    Rsites, Usites, Lsites, Dsites, onsites = [], [], [], [], []
    for r in 1:N
        nnr = [nearest_neighbours(r, L=L)...]
        push!(Rsites, nnr[1])
        push!(Usites, nnr[2])
        push!(Lsites, nnr[3])
        push!(Dsites, nnr[4])
        push!(onsites, r)
    end
    sites = [Rsites, Usites, Lsites, Dsites, onsites]

    for (b, b_sites) in enumerate(sites)
        for (d, d_sites) in enumerate(sites)
            b_idx_start = (b - 1) * N + 1
            b_idx_end = b * N
            d_idx_start = (d - 1) * N + 1
            d_idx_end = d * N

            if b == 5 || d == 5 # the on-site terms get potential V=V0
                Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V0, N=N)
            else
                Mblock = dwave_blocks(b_sites, d_sites; P=P, U=U, Uconj=Uconj, V=V1, N=N)
            end

            M[b_idx_start:b_idx_end, d_idx_start:d_idx_end] .= Mblock
        end
    end

    return M

end
