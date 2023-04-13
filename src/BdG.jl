using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations

include("../src/results.jl")

function finite_size_gap(; L::Int, t::Real, J::Real, Q::Real, μ::Real)
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ)
    E, _ = diagonalize_hamiltonian(H0)
    sort!(E)
    ΔE = [E[i+1] - E[i] for i in 1:(length(E)-1)]
end

function BdG_iteration(M::Matrix{Float64}, Δi; V0::Real, T::Real)
    M = copy(M)
    N = size(Δi)[1]
    # Fill in the off-diagonal blocks 
    M[1:N, (N+1):end] .= Δi
    M[(N+1):end, 1:N] .= conj.(Δi)
    # diagonalize this matrix to get the eigenvalues 
    E, UV = eigen(Hermitian(M))

    U = UV[1:N, :]
    V = UV[(N+1):end, :]

    # compute the gap parameter 
    @einsimd Δnew[i] := V0 / 2 * U[i, n] * V[i, n] * tanh(E[n] / (2 * T))

    return Δnew
end

function rms(a, b)
    @assert size(a) == size(b)
    norm(a - b) / prod(size(a))
end

function converge_BdG(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real, niter::Int=20, ε::Real=1e-2, tol::Union{Real,Nothing}=nothing)
    N = L * L

    # make the BdG equation matrix (fill in just diagonals)
    M = zeros(2 * N, 2 * N)
    hij = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ)
    M[1:N, 1:N] .= hij
    M[(N+1):end, (N+1):end] .= -conj.(hij)

    ## Converge Δ

    # make an initial guess for the gap parameter -- all ones along diagonal
    Δi = ε * Matrix(I, (N, N))
    Δi_diag_prev = diag(Δi)

    # iterate 
    conv = []
    for n in 1:niter
        Δi_diag = BdG_iteration(M, Δi; V0=V0, T=T) # perform one BdG iteration and get the new Δi 
        Δi = diagm(Δi_diag) # make a matrix with Δi along the diagonals 

        # calculate convergence information  
        ΔΔ = rms(Δi_diag, Δi_diag_prev) # the change in Δ from one iteration to the next 
        push!(conv, ΔΔ)
        if !isnothing(tol) && ΔΔ <= tol
            # if this change is below a certain tolerance, stop iterating and return 
            return Δi_diag, conv
        end

        # if not, keep iterating 
        Δi_diag_prev = Δi_diag
    end

    return Δi_diag_prev, conv # the converged value for Δ
end

function compute_Δ(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real, niter::Int=20, ε::Real=1e-2, tol::Union{Real,Nothing}=nothing)
    # converge the BdG 
    Δi, conv = converge_BdG(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol)

    # the gap is the maximum value of Δi 
    Δ = maximum(Δi)
end
