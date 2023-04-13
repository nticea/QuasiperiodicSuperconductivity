using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations

function expspace(start, stop, length)
    exp10.(range(start, stop, length=length))
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

function diagonalize_hamiltonian(H)
    # we must compute all eigenvalues
    vals, vecs = eigen(Hermitian(Matrix(H)))
    return vals, vecs
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

function plot_potential(; L::Int, J::Real, Q::Real)
    potmat = zeros(L, L)
    for x in 1:L
        for y in 1:L
            potmat[x, y] = aubry_andre(x, y; J=J, Q=Q)
        end
    end
    h = heatmap(potmat)
    h = xticks(h, collect(1:2:L))
    h = yticks(h, collect(1:2:L))
    h = xlabel(h, "Site (x)")
    h = ylabel(h, "Site (y)")
    h = title(h, "Potential")

    return h
end

