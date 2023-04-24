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

function noninteracting_hamiltonian(; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool=true, θ::Union{Nothing,Real}=nothing)
    # construct the kinetic part 
    Ht = square_lattice_kinetic(L=L, t=t, periodic=periodic)

    # interaction
    Hint = zeros(L * L) # this is just the diagonal 
    for x in 1:L
        for y in 1:L
            n = coordinate_to_site(x, y, L=L)
            U_xy = aubry_andre(x, y, J=J, Q=Q, L=L, θ=θ)
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

function B(; L::Int, Q::Real, θ::Real)
    # First, construct R(θ) in 2D
    c, s = cos(θ), sin(θ)
    Rθ = [c s; s -c]

    # Multiply by Q (the irrational number) to get B 
    B = Q * Rθ

    # Make periodic in L 
    BSD = round.(Int, B .* L) / L

    # Perform checks
    @assert gcd(round(Int, det(L .* BSD)), L) == 1

    return BSD
end

function aubry_andre(x, y; J::Real, Q::Real, L::Union{Int,Nothing}=nothing, θ::Union{Real,Nothing}=nothing)
    if !isnothing(L) && isnothing(θ)
        Q̃ = floor(Int, Q * L) / L
        return J * (cos(2 * π * Q̃ * (x + y)) - cos(2 * π * Q̃ * (x - y)))
    elseif isnothing(L) || isnothing(θ)
        return J * (cos(2 * π * Q * (x + y)) - cos(2 * π * Q * (x - y)))
    else
        BSD = B(L=L, Q=Q, θ=θ)
        return J * (cos(2 * π * (BSD[1, 1] * x + BSD[1, 2] * y)) + cos(2 * π * (BSD[2, 1] * x + BSD[2, 2] * y)))
    end
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

