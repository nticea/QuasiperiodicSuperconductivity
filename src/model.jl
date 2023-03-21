using LinearAlgebra
using SparseArrays
using Graphs
using ITensors
using ArnoldiMethod

function λmax(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real)
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    χ = pairfield_susceptibility(T, E=E, U=U)

    # Construct M (for s-wave, all we do is multiply χ by +V0)
    M = make_M(χ, V0)

    # Calculate Tc by finding the eigenvalues of M
    λs = diagonalize_M(M)

    return λs[1]
end

function noninteracting_hamiltonian(; L::Int, t::Real, J::Real, Q::Real, μ::Real)
    # construct the kinetic part 
    Ht = square_lattice_kinetic(L=L, t=t)

    # interaction
    Hint = zeros(L * L) # this is just the diagonal 
    for x in 1:L
        for y in 1:L
            n = coordinate_to_site(x, y, L=L)
            U_xy = aubry_andre(x, y, J=J, Q=Q)
            Hint[n] = -(U_xy + μ)
        end
    end

    Hint = spdiagm(Hint)

    return Hermitian(Ht + Hint)
end

function square_lattice_kinetic(; L::Int, t::Real, periodic::Bool=true)
    g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    H = Graphs.LinAlg.adjacency_matrix(g)
    return -t .* H
end

function coordinate_to_site(x::Int, y::Int; L::Int)
    (x - 1) * L + y
end

function aubry_andre(x, y; J::Real, Q::Real)
    J * (cos(2 * π * Q * (x + y)) - cos(2 * π * Q * (x - y)))
end

function fermi(ε::Real, T::Real)
    # note: the chemical potential was subtracted away in the Hamiltonian 
    1 / (exp(ε / T) + 1)
end

function diagonalize_hamiltonian(H)
    # we must compute all eigenvalues
    # vals, vecs = eigen(Hermitian(Matrix(H)), sortby=nothing)
    vals, vecs = eigen(Matrix(H), sortby=nothing)
    return vals, vecs
end

function diagonalize_M(M)
    decomp, _ = partialschur(M, nev=1, tol=1e-6, which=LM())
    return decomp.R
end

function pairfield_susceptibility(T; E, U)
    # χ is a Hermitian object, so we only need to compute the upper triangular 
    N = size(U)[1]
    χ = zeros(N, N, N, N)

    Threads.@threads for abcd in CartesianIndices(χ)
        χ[abcd] = χelem(Tuple(abcd)..., T=T, E=E, U=U)
    end

    return χ

end

function χelem(a::Int, b::Int, c::Int, d::Int; T::Real, E, U)
    N = size(U)[1] # this is L^2 

    χelem = 0
    for n in 1:N
        for m in 1:N
            prefac = (1 - fermi(E[n], T) - fermi(E[m], T)) / (E[n] + E[m])
            χelem += prefac * (conj(U[a, n]) * U[c, n] * conj(U[b, m]) * U[d, m] +
                               conj(U[a, n]) * U[d, n] * conj(U[b, m]) * U[c, m])
        end
    end

    return χelem
end

function make_M(χ, V0)
    M = χ * V0
    N = size(χ)[1]
    # reshape M into a 2dx2d matrix 
    a, b, c, d = Index(N, "a"), Index(N, "b"), Index(N, "c"), Index(N, "d")
    M = ITensor(M, a, b, c, d)
    # make a combiner
    ab, cd = combiner(a, b), combiner(c, d)
    M = M * ab * cd
    return array(M)
end