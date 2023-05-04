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
    H0 = Ht + Hint

    # scale the Hamiltonian
    H0 = 1 / (1 + J / 2) .* H0

    return H0
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

function nearest_neighbours(r::Int; L::Int)
    x, y = site_to_coordinate(r, L=L)

    # left 
    if x == 1
        xL = L
    else
        xL = x - 1
    end
    # right 
    if x == L
        xR = 1
    else
        xR = x + 1
    end
    # up 
    if y == L
        yU = 1
    else
        yU = y + 1
    end
    # down 
    if y == 1
        yD = L
    else
        yD = y - 1
    end

    # convert back to r representation 
    rL = coordinate_to_site(xL, y, L=L)
    rU = coordinate_to_site(x, yU, L=L)
    rR = coordinate_to_site(xR, y, L=L)
    rD = coordinate_to_site(x, yD, L=L)

    return rL, rU, rR, rD
end

function coordinate_to_site(x::Int, y::Int; L::Int)
    (x - 1) * L + y
end

function site_to_coordinate(r; L::Int)
    x = floor(Int, r / L) + 1
    y = r % L
    if y == 0
        x -= 1
        y = L
    end
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
        @error "Flipped the sign!!"
        return J * (cos(2 * π * Q̃ * (x + y)) + cos(2 * π * Q̃ * (x - y)))
    elseif isnothing(L) || isnothing(θ)
        @error "Flipped the sign!!"
        return J * (cos(2 * π * Q * (x + y)) + cos(2 * π * Q * (x - y)))
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
    xticks!(h, collect(1:2:L))
    yticks!(h, collect(1:2:L))
    xlabel!(h, "Site (x)")
    ylabel!(h, "Site (y)")
    title!(h, "Potential for J=$J")

    return h
end

