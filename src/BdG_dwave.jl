using LinearAlgebra
using SparseArrays
using Graphs
using TriangularIndices
using Einsum
using Interpolations
using Distributions

include("../src/meanfield.jl")
include("../src/model.jl")
include("../src/results.jl")

function converge_BdG_dwave(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0, Δ_init=nothing)

    # size of the BdG matrix is 2N × 2N
    N = L * L

    # make the BdG matrix (fill in just the diagonals with H0 for now)
    M = zeros(2 * N, 2 * N)
    hij = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
    M[1:N, 1:N] .= hij
    M[(N+1):end, (N+1):end] .= -conj.(hij)

    # make the interaction matrix 
    Vij = make_interaction(L=L, V0=V0, V1=V1, periodic=periodic)

    # make an initial guess for the gap parameter -- finite-size Δ gap 
    Δij_prev1 = initialize_Δ_dwave(L=L, Δ_init=Δ_init)
    Δij_prev2 = copy(Δij_prev1)

    U, V, E = nothing, nothing, nothing

    ## ITERATE ## 
    hist = [] # keep track of convergence history
    for n in 1:niter
        print(n, "-")

        # take the average of the previous 2 iterations 
        Δij_prev = (Δij_prev1 + Δij_prev2) ./ 2

        # Compute the BdG iteration 
        Δij, U, V, E = BdG_iteration_dwave(M, Δij_prev; Vij=Vij, T=T) # perform one BdG iteration and get the new Δij 

        # check for convergence 
        ΔΔ = norm(Δij .- Δij_prev)
        if ΔΔ <= tol
            Δij = ΔBdG_to_ΔLGE_flat(Δij, L=L)
            return Δij, U, V, E, hist
        end

        # keep track of convergence over time 
        push!(hist, ΔΔ)

        # add noise if desired, scaled by the largest Δij in the previous iteration (noise is set to 0 by default)
        # d = Normal(0.0, maximum(Δij) * noise)
        # ε = rand.(d, size(Δij)...)
        # Δij .+= ε

        # Update the history 
        Δij_prev2 = copy(Δij_prev1)
        Δij_prev1 = copy(Δij)
    end

    Δij = (Δij_prev1 + Δij_prev2) ./ 2
    Δij = ΔBdG_to_ΔLGE_flat(Δij, L=L)

    return Δij, U, V, E, hist # return the converged value for Δij and convergence history 
end

function compute_Δ_dwave(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0, Δ_init=nothing)

    Δij, _, _, _, hist = converge_BdG_dwave(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise)

    return Δij, hist
end


function BdG_coefficients_dwave(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real, V1::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, periodic::Bool=true, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0, Δ_init=nothing)

    # converge the BdG 
    _, U, V, E, _ = converge_BdG_dwave(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, V1=V1, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise, Δ_init=Δ_init)

    return U, V, E
end

function BdG_iteration_dwave(M::Matrix{Float64}, Δij; Vij, T::Real)
    M = copy(M)
    N = size(Δij)[1]

    # Fill in the off-diagonal blocks with Δ
    M[1:N, (N+1):end] .= Δij
    M[(N+1):end, 1:N] .= Δij' #conj.(Δij) # Check this!! 

    # diagonalize this matrix to get the eigenvalues 
    E, UV = eigen(Hermitian(M))

    # Extract the u and the v coefficients
    U = UV[1:N, :]
    V = UV[(N+1):end, :]

    # compute the updated Δij
    @einsimd Δnew[i, j] := -Vij[i, j] / 4 * (U[i, n] * conj(V[j, n]) + U[j, n] * conj(V[i, n])) * tanh(E[n] / (2 * T))

    return Δnew, U, V, E
end

function initialize_Δ_dwave(; L::Int, Δ_init=nothing)
    if isnothing(Δ_init)
        d = Uniform(0, 0.05)
        return rand.(d, L * L, L * L)
    else
        @assert size(Δ_init) == (5 * L * L,)
        return ΔLGE_to_ΔBdG(Δ_init, L=L)
    end
end

function make_interaction(; L::Int, periodic::Bool, V0::Real, V1::Real)

    # construct the adjacency matrix for the nearest-neighbours 
    g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    Vij = Graphs.LinAlg.adjacency_matrix(g)
    Vij = Matrix(Vij) # temporarily make it a dense matrix 
    Vij = convert(Matrix{Float64}, Vij)

    # the off-diagonals are multiplied by V1 
    Vij .*= V1

    # the diagonals are multiplied by V0 
    Vij[diagind(Vij)] .= V0

    return Vij
end