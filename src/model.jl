using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations

struct ModelParams
    L::Int
    t::Real
    Q::Real
    μ::Real
    θ::Real
    ϕx::Real
    ϕy::Real
    V0::Real
    V1::Real
    J::Real
    periodic::Bool
end

function ModelParams(; L, t, Q, μ, θ, ϕx, ϕy, V0, V1, J, periodic)
    ModelParams(L, t, Q, μ, θ, ϕx, ϕy, V0, V1, J, periodic)
end

function expspace(start, stop, length)
    exp10.(range(start, stop, length=length))
end

function noninteracting_hamiltonian(; L::Int, t::Real, J::Real, Q::Real, μ::Real,
    periodic::Bool=true, θ::Union{Nothing,Real}=nothing, ϕx::Real=0, ϕy::Real=0)
    # construct the kinetic part 
    Ht = square_lattice_kinetic(L=L, t=t, periodic=periodic)

    # interaction
    Hint = zeros(L * L) # this is just the diagonal 
    for x in 1:L
        for y in 1:L
            n = coordinate_to_site(x, y, L=L)
            # U_xy = aubry_andre(x, y, J=J, Q=Q, L=L, θ=θ, ϕx=ϕx, ϕy=ϕy)
            U_xy = aubry_andre(x + floor(Int, L / 2), y + floor(Int, L / 2), J=J, Q=Q, L=L, θ=θ, ϕx=ϕx, ϕy=ϕy)
            Hint[n] = -(U_xy + μ)
        end
    end
    Hint = spdiagm(Hint)
    H0 = Ht + Hint

    # scale the Hamiltonian
    #@error "I am using the J-scaled version of H"
    #H0 = 1 / (1 + J / 2) .* H0

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

function coordinate_map(; L::Int)
    mat = Matrix{Tuple{Int64,Int64}}(undef, L, L)
    for x in 1:L
        for y in 1:L
            mat[x, y] = (x, y)
        end
    end

    vec = reshape(mat, L * L)
    return mat, vec
end

function coordinate_to_site(x::Int, y::Int; L::Int)
    _, vec = coordinate_map(L=L)
    return findall(r -> r == (x, y), vec)[1]
end

function site_to_coordinate(r; L::Int)
    _, vec = coordinate_map(L=L)
    return vec[r]
end

function coordinate_to_configuration_space(x, y; L::Int, Q::Real, θ::Union{Real,Nothing})
    # prepare the vectors 
    r = [x + floor(Int, L / 2); y + floor(Int, L / 2)]
    ϕ = [ϕx; ϕy]

    # make the B matrix 
    BSD = B(L=L, Q=Q, θ=θ)
    Binv = inv(BSD)

    # compute the new ϕ
    r̃ = r + Binv * ϕ # re-parameterize r as distance between r and nearest potential min
    ϕ̃ = 2π * BSD * r̃ # find the corresponding ϕ̃
    ϕ̃ = mod.(ϕ̃, 2π)  # make ϕ periodic in 2π
end

function site_to_configuration_space(r; L::Int, Q::Real, θ::Union{Real,Nothing})
    x, y = site_to_coordinate(r, L=L)
    coordinate_to_configuration_space(x, y, L=L, Q=Q, θ=θ)
end

function B(; L::Int, Q::Real, θ::Real)
    # First, construct R(θ) in 2D
    c, s = cos(θ), sin(θ)
    Rθ = [c s; s -c]

    # Multiply by Q (the irrational number) to get B 
    Bmat = Q * Rθ

    # Make periodic in L 
    BSD = round.(Int, Bmat .* L) / L

    # Perform checks
    @assert gcd(round(Int, det(L .* BSD)), L) == 1

    return BSD
end

function aubry_andre(x, y; J::Real, Q::Real, L::Union{Int,Nothing}=nothing, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0)
    if isnothing(θ)
        # Q̃ /= √2
        Q̃ = floor(Int, Q * L) / L
        return J * (cos(2 * π * Q̃ * (x + y) + ϕx) - cos(2 * π * Q̃ * (x - y) + ϕy))
    else
        BSD = B(L=L, Q=Q, θ=θ)
        return J * (cos(2 * π * (BSD[1, 1] * x + BSD[1, 2] * y) + ϕx) + cos(2 * π * (BSD[2, 1] * x + BSD[2, 2] * y) + ϕy))
    end
    # Q̃ = floor(Int, Q * L) / L
    # return J * (cos(2 * π * Q̃ * (x + y)) - cos(2 * π * Q̃ * (x - y)))
    # if !isnothing(L) && isnothing(θ)
    #     Q̃ = floor(Int, Q * L) / L
    #     @error "Flipped the sign!!"
    #     return J * (cos(2 * π * Q̃ * (x + y)) + cos(2 * π * Q̃ * (x - y)))
    # elseif isnothing(L) || isnothing(θ)
    #     @error "Flipped the sign!!"
    #     return J * (cos(2 * π * Q * (x + y)) + cos(2 * π * Q * (x - y)))
    # else
    #     BSD = B(L=L, Q=Q, θ=θ)
    #     return J * (cos(2 * π * (BSD[1, 1] * x + BSD[1, 2] * y)) + cos(2 * π * (BSD[2, 1] * x + BSD[2, 2] * y)))
    # end
end

function fermi(ε::Real, T::Real)
    # note: the chemical potential was subtracted away in the Hamiltonian 
    1 / (exp(ε / T) + 1)
end

function plot_potential(; L::Int, J::Real, Q::Real, θ::Union{Real,Nothing}, ϕx::Real=0, ϕy::Real=0)
    potmat = zeros(L, L)
    for x in 1:L
        for y in 1:L
            potmat[x, y] = aubry_andre(x, y; L=L, J=J, Q=Q, θ=θ, ϕx=ϕx, ϕy=ϕy)
        end
    end
    h = heatmap(potmat)
    xticks!(h, collect(1:2:L))
    yticks!(h, collect(1:2:L))
    xlabel!(h, "Site (x)")
    ylabel!(h, "Site (y)")
    θ = θ_to_π(θ)
    ϕx = θ_to_π(ϕx)
    ϕy = θ_to_π(ϕy)

    title!(h, "Potential for J=$J, θ=$θ, ϕx=$ϕx, ϕy=$ϕy")

    return h
end

function finite_size_gap(; L::Int, t::Real, Q::Real, μ::Real, periodic::Bool=true, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0)
    H0 = noninteracting_hamiltonian(L=L, t=t, J=0, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
    E, _ = diagonalize_hamiltonian(H0)
    sort!(E)
    ΔE = [E[i+1] - E[i] for i in 1:(length(E)-1)]
end