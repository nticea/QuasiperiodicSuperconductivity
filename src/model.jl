using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations
import Base.copy

struct ModelParams
    L::Int
    t::Real
    Q::Real
    μ::Real
    θ::Real
    ϕx::Real
    ϕy::Real
    ϕz::Real
    V0::Real
    V1::Real
    J::Real
    periodic::Real
    ndims::Int
end

struct DiagonalizedHamiltonian
    L::Int
    t::Real
    Q::Real
    μ::Real
    θ::Real
    ϕx::Real
    ϕy::Real
    ϕz::Real
    J::Real
    periodic::Real
    ndims::Int

    # diagonalization 
    E
    U
end

function ModelParams(; L, t, Q, μ, θ, ϕx, ϕy, ϕz, V0, V1, J, periodic, ndims)
    ModelParams(L, t, Q, μ, θ, ϕx, ϕy, ϕz, V0, V1, J, periodic, ndims)
end

function copy(m::ModelParams)
    mnew = ModelParams(m.L, m.t, m.Q, m.μ, m.θ, m.ϕx, m.ϕy, m.ϕz, m.V0, m.V1, m.J, m.periodic, m.ndims)
end

function DiagonalizedHamiltonian(m::ModelParams; E, U)
    DiagonalizedHamiltonian(m.L, m.t, m.Q, m.μ, m.θ, m.ϕx, m.ϕy, m.ϕz, m.J, m.periodic, m.ndims, E, U)
end

function DiagonalizedHamiltonian(m::ModelParams)
    E, U = diagonalize_hamiltonian(m)
    DiagonalizedHamiltonian(m.L, m.t, m.Q, m.μ, m.θ, m.ϕx, m.ϕy, m.ϕz, m.J, m.periodic, m.ndims, E, U)
end

function expspace(start, stop, length)
    exp10.(range(start, stop, length=length))
end

function numsites(m::ModelParams)
    L, ndims = m.L, m.ndims
    if ndims == 2
        return L * L
    elseif ndims == 3
        return L * L * L
    else
        println("$ndims dimensions not implemented yet :((")
        return
    end
end

function noninteracting_hamiltonian(m::ModelParams; scale_model::Bool=false, shift_origin=true)
    # construct the kinetic part 
    Ht = square_lattice_kinetic(m)

    # interaction
    _, rs = coordinate_map(m)
    pot = aubry_andre.(rs, m=m, shift_origin=shift_origin)
    pot .+= μ # add chemical potential 
    Hint = spdiagm(pot)
    H0 = Ht + Hint

    # scale the Hamiltonian
    if scale_model
        @warn "I am using the J-scaled version of H"
        H0 = 1 / (1 + m.J / 2) .* H0
    end

    return H0
end

function diagonalize_hamiltonian(m; loadpath::Union{String,Nothing}=nothing)
    # we must compute all eigenvalues
    if isnothing(loadpath)
        H = noninteracting_hamiltonian(m)
        E, U = eigen(Hermitian(Matrix(H)))
        return E, U
    end

    # try to load the Hamiltonian corresponding to these parameters 
    try
        DH = load_diagonalized_H(loadpath)
        @assert DH.L == m.L && DH.t == m.t && DH.J == m.J && DH.Q == m.Q && DH.μ == m.μ && DH.θ == m.θ && DH.ϕx == m.ϕx && DH.ϕy == m.ϕy && DH.ϕz == m.ϕz && DH.periodic == m.periodic
        return DH.E, DH.U
    catch e
        @show e
        # save everything to checkpointpath 
        DH = DiagonalizedHamiltonian(m)
        save_structs(DH, loadpath)
        return DH.E, DH.U
    end
end

function square_lattice_kinetic(m::ModelParams)
    L, t, periodic, ndims = m.L, m.t, m.periodic, m.ndims

    if ndims == 2
        g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    elseif ndims == 3
        g = Graphs.SimpleGraphs.grid((L, L, L), periodic=periodic)
    else
        println("ndims=$ndims not supported (yet?)")
        return
    end

    H = Graphs.LinAlg.adjacency_matrix(g)
    return -t .* H
end

function nearest_neighbours(r::Int; m::ModelParams)
    L, ndims = m.L, m.ndims

    if ndims == 2
        x, y = site_to_coordinate(r, m=m)
    elseif ndims == 3
        x, y, z = site_to_coordinate(r, m=m)
    else
        println("$ndims dimensions not supported (yet?)")
        return
    end

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

    if ndims == 2
        # convert back to r representation 
        rL = coordinate_to_site(xL, y, m=m)
        rU = coordinate_to_site(x, yU, m=m)
        rR = coordinate_to_site(xR, y, m=m)
        rD = coordinate_to_site(x, yD, m=m)

        return rL, rU, rR, rD

    elseif ndims == 3
        # up 
        if z == L
            zU = 1
        else
            zU = z + 1
        end

        # down 
        if z == 1
            zD = L
        else
            zD = z - 1
        end

        # convert back to r representation 
        rL = coordinate_to_site(xL, y, z, m=m)
        rU = coordinate_to_site(x, yU, z, m=m)
        rR = coordinate_to_site(xR, y, z, m=m)
        rD = coordinate_to_site(x, yD, z, m=m)
        rzU = coordinate_to_site(x, y, zU, m=m)
        rzD = coordinate_to_site(x, y, zD, m=m)

        return rL, rU, rR, rD, rzU, rzD
    end
end

function coordinate_map(m::ModelParams)
    L, ndims = m.L, m.ndims

    if ndims == 2
        mat = [(x, y) for x in 1:L, y in 1:L]
        vec = reshape(mat, L * L)
    elseif ndims == 3
        mat = [(x, y, z) for x in 1:L, y in 1:L, z in 1:L]
        vec = reshape(mat, L * L * L)
    else
        println("ndims=$ndims not supported (yet?)")
        return
    end

    return mat, vec
end

function coordinate_to_site(x::Int, y::Int; m::ModelParams)
    @assert m.ndims == 2
    _, vec = coordinate_map(m)
    return findall(r -> r == (x, y), vec)[1]
end

function coordinate_to_site(x::Int, y::Int, z::Int; m::ModelParams)
    @assert m.ndims == 3
    _, vec = coordinate_map(m)
    return findall(r -> r == (x, y, z), vec)[1]
end

function site_to_coordinate(r; m::ModelParams)
    _, vec = coordinate_map(m)
    return vec[r]
end

function site_to_configuration_space(r0; m::ModelParams, shift_origin::Bool=true)
    L, ndims = m.L, m.ndims
    ϕx, ϕy, ϕz = m.ϕx, m.ϕy, m.ϕz

    # prepare the vectors 
    if ndims == 2
        x, y = site_to_coordinate(r0, m=m)

        if shift_origin
            x += floor(Int, L / 2)
            y += floor(Int, L / 2)
        end

        r = [x; y]
        ϕ = [ϕx; ϕy]

    elseif ndims == 3
        x, y, z = site_to_coordinate(r0, m=m)

        if shift_origin
            x += floor(Int, L / 2)
            y += floor(Int, L / 2)
            z += floor(Int, L / 2)
        end

        r = [x; y; z]
        ϕ = [ϕx; ϕy; ϕz]

    else
        println("$ndims dimensions not supported")
    end

    # make the B matrix 
    BSD = Bmatrix(m)
    Binv = inv(BSD)

    # compute the new ϕ
    r̃ = r + Binv * ϕ # re-parameterize r as distance between r and nearest potential min
    ϕ̃ = 2π * BSD * r̃ # find the corresponding ϕ̃
    ϕ̃ = mod.(ϕ̃, 2π)  # make ϕ periodic in 2π
end

function Bmatrix(m::ModelParams)
    L, Q, θ, ndims = m.L, m.Q, m.θ, m.ndims
    c, s = cos(θ), sin(θ)

    #Rθ = [c s; s -c]
    if ndims == 2
        Rθ = [c -s; s c]
    elseif ndims == 3
        Rθ = [c -s 0; s c 0; 0 0 1] # this is around y axis 
    # Rθ = [1 0 0; 0 c -s; 0 s c] this is around x axis 
    else
        println("ndims=$ndims not supported (yet?)")
        return
    end

    # Multiply by Q (the irrational number) to get B 
    Bmat = Q * Rθ

    # Make periodic in L 
    BSD = round.(Int, Bmat .* L) / L

    # Perform checks
    @assert gcd(round(Int, det(L .* BSD)), L) == 1

    return BSD
end

function aubry_andre(xy; m::ModelParams, shift_origin::Bool=true)
    J, Q, L, θ, ϕx, ϕy, ϕz, ndims = m.J, m.Q, m.L, m.θ, m.ϕx, m.ϕy, m.ϕz, m.ndims

    if shift_origin
        xy = [a + floor(Int, L / 2) for a in xy]
    end

    # make the rotation matrix 
    BSD = Bmatrix(m)

    if ndims == 2
        x, y = xy
        t1 = 2 * π * (BSD[1, 1] * x + BSD[1, 2] * y) + ϕx
        t2 = 2 * π * (BSD[2, 1] * x + BSD[2, 2] * y) + ϕy
        return J * (cos(t1) + cos(t2))
    elseif ndims == 3
        x, y, z = xy
        t1 = 2 * π * (BSD[1, 1] * x + BSD[1, 2] * y + BSD[1, 3] * z) + ϕx
        t2 = 2 * π * (BSD[2, 1] * x + BSD[2, 2] * y + BSD[2, 3] * z) + ϕy
        t3 = 2 * π * (BSD[3, 1] * x + BSD[3, 2] * y + BSD[3, 3] * z) + ϕz
        return J * (cos(t1) + cos(t2) + cos(t3))
    else
        println("ndims=$ndims not supported (yet?)")
        return
    end
end

function fermi(ε::Real, T::Real)
    # note: the chemical potential was subtracted away in the Hamiltonian 
    1 / (exp(ε / T) + 1)
end

function plot_potential(m::ModelParams; slice::Int=1)
    L, J, θ, ϕx, ϕy, μ = m.L, m.J, m.θ, m.ϕx, m.ϕy, m.μ

    _, rs = coordinate_map(m)
    potmat = aubry_andre.(rs, m=m)
    potmat .+= μ # add chemical potential
    if ndims == 2
        potmat = reshape(potmat, L, L)
    elseif ndims == 3
        potmat = reshape(potmat, L, L, L)
        potmat = potmat[:, :, slice]
    end

    numpts = 10
    cm = cgrad(:bwr, 2 * numpts + 1, categorical=true)

    function colour_gradient(x1::Int, x2::Int; arr)
        val = arr[x1, x2]
        max = 3 * J
        idx = floor(Int, val / max * numpts + numpts + 1)
        return cm[idx]
    end

    plot(xaxis=" ", yaxis=" ", legend=false)
    heatmap([-2*J 2; 3 2*J], color=cm, visible=false)  # Dummy heatmap to generate colorbar

    # Create the colorbar

    h = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
    for x in 1:L
        for y in 1:L
            # onsite dot 
            scatter!(h, [x], [y], ms=10, c=colour_gradient(x, y, arr=potmat), legend=:false, aspect_ratio=:equal)
        end
    end

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

function momentum_components(m::ModelParams; r::Int)
    ndims = m.ndims
    if ndims == 2
        x, y = site_to_coordinate(r, m=m)
        return (2 * π * x / L, 2 * π * y / L)
    elseif ndims == 3
        x, y, z = site_to_coordinate(r, m=m)
        return (2 * π * x / L, 2 * π * y / L, 2 * π * z / L)
    else
        println("$dims dimensions not implemented")
    end
end

function finite_size_gap(m::ModelParams)
    H0 = noninteracting_hamiltonian(m)
    E, _ = diagonalize_hamiltonian(H0)
    sort!(E)
    ΔE = [E[i+1] - E[i] for i in 1:(length(E)-1)]
end