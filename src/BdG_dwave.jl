using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations
using ITensors

include("../src/results.jl")

function finite_size_gap(; L::Int, t::Real, Q::Real, μ::Real, periodic::Bool)
    H0 = noninteracting_hamiltonian(L=L, t=t, J=0, Q=Q, μ=μ, periodic=periodic)
    E, _ = diagonalize_hamiltonian(H0)
    sort!(E)
    ΔE = [E[i+1] - E[i] for i in 1:(length(E)-1)]
end

function BdG_iteration(M::Matrix{Float64}, Δi; Vij, T::Real)
    M = copy(M)
    N = size(Δi)[1]
    # Fill in the off-diagonal blocks 
    M[1:N, (N+1):end] .= Δi
    M[(N+1):end, 1:N] .= conj.(Δi)

    # diagonalize this matrix to get the eigenvalues 
    E, UV = eigen(Hermitian(M))

    # check this !!!
    U = UV[1:N, :]
    V = UV[(N+1):end, :]

    # compute the gap parameter 
    @einsimd Δnew[i, j] := Vij[i, j] / 2 * U[i, n] * conj(V[j, n]) * tanh(E[n] / (2 * T))

    # # try this using ITensors
    # # define the indices
    # dim_i, dim_j = size(Vij)
    # dim_n = length(E)
    # i, j, n = Index(dim_i, "i"), Index(dim_j, "j"), Index(dim_n, "n")

    # # make the tensors 
    # Vij = ITensor(Matrix(Vij), i'', j'')
    # U = ITensor(U, i', n)
    # conjV = ITensor(conj.(V), j', n)
    # tanhE = tanh.(E) ./ (2 * T)
    # tanhE = ITensor(tanhE, n)

    # # perform the contraction 
    # println("New way...")
    # @time Δnew = 1 / 2 * Vij * ITensors.δ(i', i'', i) * U * ITensors.δ(j', j'', j) #* conjV * tanhE
    # @show inds(Δnew)
    # Δnew = array(Δnew)

    return Δnew
end

function rms(a, b)
    @assert size(a) == size(b)
    norm(a - b) / prod(size(a))
end

function initialize_Δ(; L::Int, t::Real, Q::Real, μ::Real, V0::Real, V1::Real, periodic::Bool=true)

    # make an initial guess for the gap parameter -- size of fs gap 
    fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, periodic=periodic))
    @show fsgap

    if V1 == 0
        return sparse(diagm(fsgap * ones(L * L)))
    end

    g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    adjm = Graphs.LinAlg.adjacency_matrix(g)
    adjm = Matrix(adjm) # temporarily make it a dense matrix 
    adjm = convert(Matrix{Float64}, adjm)

    Δij = fsgap .* adjm

    return sparse(Δij)
end

function make_interaction(; L::Int, periodic::Bool, V0::Real, V1::Real)
    g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    Vij = Graphs.LinAlg.adjacency_matrix(g)

    Vij = Matrix(Vij) # temporarily make it a dense matrix 
    Vij = convert(Matrix{Float64}, Vij)

    # the off-diagonals are multiplied by V1 
    Vij .*= V1

    # the diagonals are multiplied by V0 
    Vij[diagind(Vij)] .= V0

    return sparse(Vij)
end

function converge_BdG(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, niter::Int=100, tol::Union{Real,Nothing}=nothing, θ::Union{Real,Nothing})
    N = L * L

    # make the BdG equation matrix (fill in just diagonals)
    M = zeros(2 * N, 2 * N)
    hij = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, periodic=periodic)
    M[1:N, 1:N] .= hij
    M[(N+1):end, (N+1):end] .= -conj.(hij)

    # make the Vij matrix 
    Vij = make_interaction(L=L, V0=V0, V1=V1, periodic=periodic)

    # make an initial guess for the gap parameter -- finite-size Δ gap  
    Δij_prev = initialize_Δ(L=L, t=t, Q=Q, μ=μ, periodic=periodic, V0=V0, V1=V1)

    # iterate 
    conv = []
    max_Δ = []
    for n in 1:niter
        Δij = BdG_iteration(M, Δij_prev; Vij=Vij, T=T) # perform one BdG iteration and get the new Δi 
        @show maximum(Δij)

        # calculate convergence information  
        ΔΔ = abs(maximum(abs.(Δij)) - maximum(abs.(Δij_prev)))
        push!(conv, (ΔΔ / maximum(Δij)))
        push!(max_Δ, maximum(Δij))
        @show conv[end]
        # if this change is below a certain tolerance, stop iterating and return 
        if !isnothing(tol) && (ΔΔ / maximum(Δij) <= tol) && (n > 1)
            return Δij, conv, max_Δ
        end

        # if not, keep iterating 
        Δij_prev = Δij
    end

    return Δij_prev, conv, max_Δ # the converged value for Δ
end

function compute_Δ_dwave(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real, V1::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, θ::Union{Real,Nothing}, periodic::Bool)
    # converge the BdG 
    Δij, conv, max_Δ = converge_BdG(T, L=L, t=t, J=J, Q=Q, μ=μ, periodic=periodic, V0=V0, V1=V1, tol=tol, θ=θ, niter=niter)
end
