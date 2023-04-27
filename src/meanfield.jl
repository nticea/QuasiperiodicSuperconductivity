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

function λmax(T; L::Int, t::Real, J::Real, Q::Real, θ::Union{Real,Nothing}=nothing, μ::Real, V0::Real, V1::Real=0, symmetry::String="s-wave")
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    if symmetry == "s-wave" || symmetry == "d"
        M = swave(T, E=E, U=U, V0=V0)
    elseif symmetry == "d-wave" || symmetry == "d"
        M = dwave(T, L=L, E=E, U=U, V0=V0, V1=V1)
    else
        @error "Pairing symmetry not recognized"
        return
    end

    # Calculate Tc by finding the eigenvalues of M
    λs = diagonalize_M(M)

    ## TODO: DIAGONALIZE EACH BLOCK INDEPENDENTLY

    return M, λs[1]
end

function diagonalize_M(M)
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

function dwave(T::Real; L, E, U, V0, V1)
    println("d-wave configuration")
    N = L^2
    Vs = [V1, V1, V1, V1, V0]

    # make the prefactor
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    P = (1 .- fnm) ./ Enm

    # because a=b always, we can do this multiplication quickly 
    Uconj = conj.(U)
    @einsimd PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]

    ## TODO:: CONSTRUCT EACH BLOCK INDEPENDENTLY

    # put together the M matrix 
    M = zeros(N, 5, N, 5)
    Threads.@threads for r in 1:N
        for rp in 1:N
            nn = [nearest_neighbours(rp, L=L)...]
            push!(nn, rp)
            for (idx, rpd) in enumerate(nn)
                PUUr = PUU[r, :, :]
                Ucn = U[rp, :]
                Udm = U[rpd, :]
                @einsimd UU[m, n] := Ucn[n] * Udm[m]

                @einsimd χ := PUUr[n, m] * UU[m, n]

                M[r, idx, rp, idx] = Vs[idx] * χ
            end
        end
    end

    # reshape M 
    M = reshape(M, 5 * N, 5 * N)

    return M
end
