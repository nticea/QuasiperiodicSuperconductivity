using LinearAlgebra, StatsBase
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations
import Base.copy
using FFTW

include("utilities.jl")

struct ModelParams
    L::Int
    t::Real
    Q::Real
    μ::Real
    θ::Real
    ϕx::Real
    ϕy::Real
    ϕz::Real
    ϕBx::Real
    ϕBy::Real
    ϕBz::Real
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
    ϕBx::Real
    ϕBy::Real
    ϕBz::Real
    J::Real
    periodic::Real
    ndims::Int

    # diagonalization 
    E
    U
end

function ModelParams(; L, t, Q, μ, θ, ϕx, ϕy, ϕz, V0, V1, J, periodic, ndims, ϕBx=0, ϕBy=0, ϕBz=0)
    ModelParams(L, t, Q, μ, θ, ϕx, ϕy, ϕz, ϕBx, ϕBy, ϕBz, V0, V1, J, periodic, ndims)
end

function copy(m::ModelParams)
    mnew = ModelParams(m.L, m.t, m.Q, m.μ, m.θ, m.ϕx, m.ϕy, m.ϕz, m.ϕBx, m.ϕBy, m.ϕBz, m.V0, m.V1, m.J, m.periodic, m.ndims)
end

function DiagonalizedHamiltonian(m::ModelParams; E, U)
    DiagonalizedHamiltonian(m.L, m.t, m.Q, m.μ, m.θ, m.ϕx, m.ϕy, m.ϕz, m.ϕBx, m.ϕBy, m.ϕBz, m.J, m.periodic, m.ndims, E, U)
end

function DiagonalizedHamiltonian(m::ModelParams; loadpath::Union{String,Nothing}=nothing)
    E, U = diagonalize_hamiltonian(m, loadpath=loadpath)
    DiagonalizedHamiltonian(m.L, m.t, m.Q, m.μ, m.θ, m.ϕx, m.ϕy, m.ϕz, m.ϕBx, m.ϕBy, m.ϕBz, m.J, m.periodic, m.ndims, E, U)
end

function expspace(start, stop, length)
    exp10.(range(start, stop, length=length))
end

function numsites(m::Union{ModelParams,DiagonalizedHamiltonian})
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
    if m.ϕBx == 0 && m.ϕBy == 0 && m.ϕBz == 0
        Ht = square_lattice_kinetic(m)
    else
        Ht = lattice_with_flux(m, ϕ1=m.ϕBx, ϕ2=m.ϕBy, ϕ3=m.ϕBz)
    end

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
        @assert DH.L == m.L && DH.t == m.t && DH.J == m.J && DH.Q == m.Q && DH.μ == m.μ && DH.θ == m.θ && DH.ϕx == m.ϕx && DH.ϕy == m.ϕy && DH.ϕz == m.ϕz && DH.ϕBx == m.ϕBx && DH.ϕBy == m.ϕBy && DH.ϕBz == m.ϕBz && DH.periodic == m.periodic && DH.ndims == m.ndims
        println("Loading a pre-diagonalized Hamiltonian")
        return DH.E, DH.U
    catch e
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

"""
NOTE: ϕᵢ are the magnetic fluxes being threaded through each unit cell,
NOT the potential offsets 
"""
function lattice_with_flux(m::ModelParams; ϕ1::Real=0, ϕ2::Real=0, ϕ3::Real=0)
    ndims, t = m.ndims, m.t
    nsites = numsites(m)
    H = zeros(ComplexF64, nsites, nsites)  # Initialize Hamiltonian matrix
    nni = nearest_neighbours.(collect(1:nsites), m=m)
    if ndims == 3
        for i in 1:nsites
            rL, rU, rR, rD, rzU, rzD = nni[i]
            for j in 1:nsites
                if j == rL || j == rR
                    H[i, j] = -t * exp(1im * ϕ1)
                elseif j == rU || j == rD
                    H[i, j] = -t * exp(1im * ϕ2)
                elseif j == rzU || j == rzD
                    H[i, j] = -t * exp(1im * ϕ3)
                end
            end
        end
    elseif ndims == 2
        for i in 1:nsites
            rL, rU, rR, rD = nni[i]
            for j in 1:nsites
                rL, rU, rR, rD, = nearest_neighbours(i, m=m)
                if j == rL || j == rR
                    H[i, j] = -t * exp(1im * ϕ1)
                elseif j == rU || j == rD
                    H[i, j] = -t * exp(1im * ϕ2)
                end
            end
        end
    else
        println("the worst")
    end
    return sparse(H)
end

function nearest_neighbours(r::Int; m::ModelParams)
    L, ndims, periodic = m.L, m.ndims, m.periodic

    if ndims == 2
        x, y = site_to_coordinate(r, m=m)
    elseif ndims == 3
        x, y, z = site_to_coordinate(r, m=m)
    else
        println("$ndims dimensions not supported (yet?)")
        return
    end

    # left 
    if x == 1 && periodic
        xL = L
    else
        xL = x - 1
    end
    # right 
    if x == L && periodic
        xR = 1
    else
        xR = x + 1
    end
    # up 
    if y == L && periodic
        yU = 1
    else
        yU = y + 1
    end
    # down 
    if y == 1 && periodic
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
        if z == L && periodic
            zU = 1
        else
            zU = z + 1
        end

        # down 
        if z == 1 && periodic
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
    idxs = findall(r -> r == (x, y), vec)
    if length(idxs) > 0
        return idxs[1]
    else
        return nothing
    end
end

function coordinate_to_site(x::Int, y::Int, z::Int; m::ModelParams)
    @assert m.ndims == 3
    _, vec = coordinate_map(m)
    idxs = findall(r -> r == (x, y, z), vec)
    if length(idxs) > 0
        return idxs[1]
    else
        return nothing
    end
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
        # Rθ = [c -s 0; s c 0; 0 0 1] # this is around y axis 
        # Rθ = [1 0 0; 0 c -s; 0 s c] this is around x axis 
        # Rθ = [c^2+s^3 c*s c*s^2-c*s
        #     c*s -s c^s
        #     c*s^2-c*s c^2 c^2*s+s^2]
        Rθ = [1 1 0; 0 1 1; 1 0 1]
        #println("Rotation matrix with no C₄ symmetry!")
    else
        println("ndims=$ndims not supported (yet?)")
        return
    end

    # Multiply by Q (the irrational number) to get B 
    Bmat = Q * Rθ

    # Make periodic in L
    if m.periodic
        BSD = round.(Int, Bmat .* L) / L
        # Perform checks
        @assert gcd(round(Int, det(L .* BSD)), L) == 1
    else
        BSD = Bmat
    end

    return BSD
end

function aubry_andre(xy; m::ModelParams, shift_origin::Bool=true, normalize_SD_to_1::Bool=true)
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

function density_of_states(H::DiagonalizedHamiltonian; nbins::Int)
    return hist_counts(H.E, nbins=nbins, normalize=true)
    # h = fit(Histogram, E, nbins=nbins)
    # return h.weights
end

function IPR_momentum(H::DiagonalizedHamiltonian)
    # transform to momentum space
    Uq = fourier_transform_eigenstates(H)

    # IPR = ∑ₖ |ψ(k)|⁴ 
    function _IPR_calc(u)
        u2 = u .* conj.(u)
        u4 = u2 .^ 2
        return real(sum(u4))
    end

    ipr = _IPR_calc.(eachcol(Uq))
    return ipr
end

function IPR_real(H::DiagonalizedHamiltonian)
    U = H.U

    # IPR = ∑ₓ |ψ(x)|⁴ 
    function _IPR_calc(u)
        u2 = u .* conj.(u)
        u4 = u2 .^ 2
        return real(sum(u4))
    end

    ipr = _IPR_calc.(eachcol(U))
    return ipr
end

function IPR_momentum(m::ModelParams)
    DH = DiagonalizedHamiltonian(m)
    return IPR_momentum(DH)
end

function IPR_real(m::ModelParams)
    DH = DiagonalizedHamiltonian(m)
    return IPR_real(DH)
end

function mean_level_spacing(H::DiagonalizedHamiltonian)
    E = H.E

    # get spacings
    δns = zeros(length(E) - 1)
    for i in 2:length(E)
        δns[i-1] = E[i] - E[i-1]
    end

    # take ratios 
    rs = []
    for i in 2:length(δns)
        r = minimum([δns[i], δns[i-1]]) / maximum([δns[i], δns[i-1]])
        push!(rs, r)
    end

    return rs
end

function mean_level_spacing(m::ModelParams)
    DH = DiagonalizedHamiltonian(m)
    return mean_level_spacing(DH)
end

function closest_indices(arr, value, n)
    # Calculate the absolute differences between each element and the value
    abs_diff = abs.(arr .- value)
    # Sort the indices of the array based on the absolute differences
    sorted_indices = sortperm(abs_diff)
    # Take the first `n` indices from the sorted list
    return sorted_indices[1:n]
end

function multifractal_mean(m::ModelParams; E₀::Real, λ::Real, loadpath::Union{String,Nothing}=nothing, num_avg::Int=10)
    DH = DiagonalizedHamiltonian(m; loadpath=loadpath)
    return multifractal_mean(DH, E₀=E₀, λ=λ, num_avg=num_avg)
end

function multifractal_mean(m::DiagonalizedHamiltonian; E₀::Real, λ::Real, num_avg::Int=10)
    ndims, L = m.ndims, m.L

    ℓ = Int(λ * L)
    E, U = m.E, m.U
    sortidx = sortperm(E)
    E = E[sortidx]
    U = U[:, sortidx]

    idxs = closest_indices(E, E₀, num_avg)
    us = U[:, idxs]

    # iterate this for each column of U 
    function get_α̃(u)
        function cg_weight(uᵢ)
            return sum(uᵢ .* conj.(uᵢ))
        end

        if ndims == 2
            u = reshape(u, L, L)
            ucubes = split_into_squares(u, ℓ)
        elseif ndims == 3
            u = reshape(u, L, L, L)
            ucubes = split_into_cubes(u, ℓ)
        end

        # then compute the sum within each square 
        uvec = reshape(ucubes, prod(size(ucubes)))
        μk = [cg_weight(uᵢ) for uᵢ in uvec]
        α̃ = log.(μk) ./ log(λ)
        return mean(α̃)
    end

    all_α̃s = get_α̃.(eachcol(us))
    return mean(all_α̃s)
end

# We need to transform each of the eigenvectors into 2D space! 
function fourier_transform_eigenstates(DH; minus=false)

    ndims, L, U = DH.ndims, DH.L, DH.U

    function FT_vec(u)
        # first, map back to 2D space 
        if ndims == 2
            N = L * L
            u = reshape(u, L, L)
        elseif ndims == 3
            N = L * L * L
            u = reshape(u, L, L, L)
        else
            println("$dims dimensions not implemented")
            return
        end

        # take FT along all spatial dimensions 
        if !minus # this is the expression for U_{q}
            uq = fft(u)
        else # this is the expression for U_{-q}
            uq = conj.(fft(conj.(u)))
        end

        # reshape it back and normalize. FFTW does not normalize!!
        uq = reshape(uq, N) ./ √N

        return uq
    end

    Uq = FT_vec.(eachcol(U))
    return hcat(Uq...)
end

function compute_scaling_properties(m::ModelParams; λ::Real, loadpath::Union{String,Nothing}=nothing, num_avg::Int=10)
    H = DiagonalizedHamiltonian(m, loadpath=loadpath)

    # multifractal mean
    α₀ = multifractal_mean(H; E₀=E₀, λ=λ, num_avg=num_avg)

    # IPR 
    idxs = closest_indices(H.E, E₀, num_avg)
    # calculate the real IPR 
    ipr = IPR_real(H)
    ipr_real = mean(ipr[idxs])

    # calculate the real IPR 
    ipr = IPR_momentum(H)
    ipr_k = mean(ipr[idxs])

    return α₀, ipr_real, ipr_k
end

