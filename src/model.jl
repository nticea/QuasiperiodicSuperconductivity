using LinearAlgebra
using SparseArrays
using Graphs
using ITensors
using ArnoldiMethod

function λmax(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real)
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ)
    #heatmap(Matrix(H0),yflip=true,clims=(-maximum(abs.(H0)),maximum(abs.(H0))),cmap=:bwr)

    # Diagonalize this Hamiltonian
    E, U = @time diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    χ = @time pairfield_singlet(T, E=E, U=U)

    # Construct M (for s-wave, all we do is multiply χ by +V0)
    M = make_M(χ, V0)

    # Calculate Tc by finding the eigenvalues of M
    λs = @time diagonalize_M(M)

    return λs[1]
end

function noninteracting_hamiltonian(; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool=true)
    # construct the kinetic part 
    Ht = square_lattice_kinetic(L=L, t=t, periodic=periodic)

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

    return Ht + Hint
end

function square_lattice_kinetic(; L::Int, t::Real, periodic::Bool=true)
    g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    H = Graphs.LinAlg.adjacency_matrix(g)
    return -t .* H
end

function coordinate_to_site(x::Int, y::Int; L::Int)
    (x - 1) * L + y
end

function site_to_coordinate(r; L::Int)
    x = floor(Int, r / L) + 1
    y = r % L
    return x, y
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
    vals, vecs = eigen(Hermitian(Matrix(H)))
    return vals, vecs
end

function diagonalize_M(M)
    decomp, _ = partialschur(M, nev=1, tol=1e-6, which=LM())
    return decomp.R
end

function pairfield_singlet(T::Real; E, U)
    N = size(U)[1]
    χ = zeros(N, N)

    # make the prefactor
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    Pnm = (1 .- fnm) ./ Enm

    # upper triangular part 
    for r in 1:N
        for rprime in (r+1):N
            χ[r, rprime] = χelem(r, r, rprime, rprime; U=U, P=Pnm)
        end
    end

    # lower triangular part
    χ += χ'

    # diagonal 
    for r in 1:N
        χ[r, r] = χelem(r, r, r, r; U=U, P=Pnm)
    end

    return χ
end

function χelem(a::Int, b::Int, c::Int, d::Int; U, P)
    N = size(U)[1] # this is L^2 

    χelem = 0
    for n in 1:N
        for m in 1:N
            χelem += P[n, m] * (conj(U[a, n]) * U[c, n] * conj(U[b, m]) * U[d, m] +
                                conj(U[a, n]) * U[d, n] * conj(U[b, m]) * U[c, m])
        end
    end

    return 1 / 2 * χelem
end

function make_M(χ, V0)
    return χ * V0 # scale by V0
end