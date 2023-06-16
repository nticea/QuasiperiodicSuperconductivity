# using LinearAlgebra
# using SparseArrays
# using Graphs
# using TriangularIndices
# using Einsum
# using Interpolations
# using Distributions

# include("../src/meanfield.jl")
# include("../src/model.jl")

# function compute_Δ_dwave(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, θ::Union{Real,Nothing},
#     ϕx::Real=0, ϕy::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0)

#     # size of the BdG matrix is 2N × 2N
#     N = L * L

#     # make the BdG matrix (fill in just the diagonals with H0 for now)
#     M = zeros(2 * N, 2 * N)
#     hij = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
#     M[1:N, 1:N] .= hij
#     M[(N+1):end, (N+1):end] .= -conj.(hij)

#     # make the interaction matrix 
#     Vij = make_interaction(L=L, V0=V0, V1=V1, periodic=periodic)

#     # make an initial guess for the gap parameter -- finite-size Δ gap 
#     Δij_prev1 = initialize_Δ_dwave(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
#     Δij_prev2 = copy(Δij_prev1)

#     ## ITERATE ## 
#     max_Δ = [] # keep track of convergence history 
#     for n in 1:niter
#         print(n, "-")

#         # the Δ_ij to use is the average of the previous two Interpolations
#         Δij_prev = Δij_prev1#(Δij_prev1 + Δij_prev2) ./ 2

#         # Perform the BdG iteration
#         Δij = BdG_iteration_dwave(M, Δij_prev; Vij=Vij, T=T) # perform one BdG iteration and get the new Δij 

#         ΔΔ = norm(Δij .- Δij_prev) ./ (2N * 2N)
#         if ΔΔ < tol
#             return Δij_prev, max_Δ
#         end

#         # add noise if desired, scaled by the largest Δij in the previous iteration (noise is set to 0 by default)
#         if n < 20
#             d = Normal(0.0, maximum(Δij) * noise)
#             ε = rand.(d, size(Δij)...)
#             Δij .+= ε
#         end

#         # Update the history 
#         Δij_prev2 = copy(Δij_prev1)
#         Δij_prev1 = copy(Δij)

#         # keep track of convergence over time 
#         push!(max_Δ, maximum(Δij))
#     end

#     Δij_prev = (Δij_prev1 + Δij_prev2) / 2

#     return Δij_prev, max_Δ # return the converged value for Δij and convergence history 
# end

# function BdG_iteration_dwave(M::Matrix{Float64}, Δij; Vij, T::Real)
#     M = copy(M)
#     N = size(Δij)[1]

#     # Fill in the off-diagonal blocks with Δ
#     M[1:N, (N+1):end] .= Δij
#     M[(N+1):end, 1:N] .= Δij' #conj.(Δij) # Check this!! 

#     # diagonalize this matrix to get the eigenvalues 
#     E, UV = eigen(Hermitian(M))

#     # Extract the u and the v coefficients
#     U = UV[1:N, :]
#     V = UV[(N+1):end, :]

#     # compute the updated Δij
#     @einsimd Δnew[i, j] := -Vij[i, j] / 4 * (U[i, n] * conj(V[j, n]) + U[j, n] * conj(V[i, n])) * tanh(E[n] / (2 * T))

#     return Δnew
# end

# function initialize_Δ_dwave(; L::Int, t::Real, Q::Real, μ::Real, θ::Union{Real,Nothing}, ϕx::Real, ϕy::Real, periodic::Bool=true)

#     # make an initial guess for the gap parameter -- size of fs gap 
#     fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic))

#     # construct the adjacency matrix for the nearest-neighbours 
#     g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
#     Δij = Graphs.LinAlg.adjacency_matrix(g)
#     Δij = Matrix(Δij) # temporarily make it a dense matrix 
#     Δij = convert(Matrix{Float64}, Δij)

#     # the off-diagonals (nearest-neighbours) are fsgap size 
#     Δij .*= fsgap

#     # the diagonals (on-site) are zero
#     # Δij[diagind(Δij)] .= fsgap

#     return Δij
# end

# function make_interaction(; L::Int, periodic::Bool, V0::Real, V1::Real)

#     # construct the adjacency matrix for the nearest-neighbours 
#     g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
#     Vij = Graphs.LinAlg.adjacency_matrix(g)
#     Vij = Matrix(Vij) # temporarily make it a dense matrix 
#     Vij = convert(Matrix{Float64}, Vij)

#     # the off-diagonals are multiplied by V1 
#     Vij .*= V1

#     # the diagonals are multiplied by V0 
#     Vij[diagind(Vij)] .= V0

#     return Vij
# end

# function finite_size_gap(; L::Int, t::Real, Q::Real, μ::Real, θ::Union{Real,Nothing}, ϕx::Real, ϕy::Real, periodic::Bool)
#     H0 = noninteracting_hamiltonian(L=L, t=t, J=0, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
#     E, _ = diagonalize_hamiltonian(H0)
#     sort!(E)
#     ΔE = [E[i+1] - E[i] for i in 1:(length(E)-1)]
# end

using LinearAlgebra
using SparseArrays
using Graphs
using TriangularIndices
using Einsum
using Interpolations
using Distributions

include("../src/meanfield.jl")
include("../src/model.jl")

function compute_Δ_dwave(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0)

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
    Δij_prev = initialize_Δ_dwave(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)

    ## ITERATE ## 
    max_Δ = [] # keep track of convergence history 
    for n in 1:niter
        print(n, "-")
        Δij = BdG_iteration_dwave(M, Δij_prev; Vij=Vij, T=T) # perform one BdG iteration and get the new Δij 

        # check for convergence 
        ΔΔ = norm(Δij .- Δij_prev)
        if ΔΔ <= tol
            return Δij_prev, max_Δ
        end

        # add noise if desired, scaled by the largest Δij in the previous iteration (noise is set to 0 by default)
        d = Normal(0.0, maximum(Δij) * noise)
        ε = rand.(d, size(Δij)...)
        Δij .+= ε

        # keep track of convergence over time 
        push!(max_Δ, maximum(Δij))

        Δij_prev = copy(Δij)
    end

    return Δij_prev, max_Δ # return the converged value for Δij and convergence history 
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

    return Δnew
end

function initialize_Δ_dwave(; L::Int, t::Real, Q::Real, μ::Real, θ::Union{Real,Nothing}, ϕx::Real, ϕy::Real, periodic::Bool=true)

    # # make an initial guess for the gap parameter -- size of fs gap 
    # fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic))

    # # construct the adjacency matrix for the nearest-neighbours 
    # g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    # Δij = Graphs.LinAlg.adjacency_matrix(g)
    # Δij = Matrix(Δij) # temporarily make it a dense matrix 
    # Δij = convert(Matrix{Float64}, Δij)

    # # the off-diagonals (nearest-neighbours) are fsgap size 
    # Δij .*= 0.05#fsgap

    # # the diagonals (on-site) are zero
    # #Δij[diagind(Δij)] .= fsgap

    # return Δij
    N = L * L
    d = Uniform(0, 0.05)
    Δij = rand.(d, N, N)
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

function finite_size_gap(; L::Int, t::Real, Q::Real, μ::Real, θ::Union{Real,Nothing}, ϕx::Real, ϕy::Real, periodic::Bool)
    H0 = noninteracting_hamiltonian(L=L, t=t, J=0, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
    E, _ = diagonalize_hamiltonian(H0)
    sort!(E)
    ΔE = [E[i+1] - E[i] for i in 1:(length(E)-1)]
end
