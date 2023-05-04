using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations

include("../src/results.jl")

function finite_size_gap(; L::Int, t::Real, Q::Real, μ::Real, periodic::Bool=true)
    H0 = noninteracting_hamiltonian(L=L, t=t, J=0, Q=Q, μ=μ, periodic=periodic)
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

    # check this !!!
    U = UV[1:N, :]
    V = UV[(N+1):end, :]

    @einsimd Δnew[i] := V0 / 2 * U[i, n] * conj(V[i, n]) * tanh(E[n] / (2 * T))

    return Δnew
end

function rms(a, b)
    @assert size(a) == size(b)
    norm(a - b) / prod(size(a))
end

function converge_BdG(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool,
    V0::Real, niter::Int=100, tol::Union{Real,Nothing}=nothing,
    θ::Union{Real,Nothing}, fsgap::Union{Real,Nothing}=nothing)

    N = L * L

    # make the BdG equation matrix (fill in just block diagonals)
    M = zeros(2 * N, 2 * N)
    hij = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, periodic=periodic)
    M[1:N, 1:N] .= hij
    M[(N+1):end, (N+1):end] .= -conj.(hij)

    if isnothing(fsgap)
        fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, periodic=periodic))
    end

    # initial gues for Δ_i
    Δi = fsgap * Matrix(I, (N, N))
    Δi_diag_prev = diag(Δi)

    # iterate 
    conv = []
    max_Δ = []
    for n in 1:niter
        Δi_diag = BdG_iteration(M, Δi; V0=V0, T=T) # perform one BdG iteration and get the new Δi 
        Δi = diagm(Δi_diag) # make a matrix with Δi along the diagonals

        # calculate convergence information  
        ΔΔ = abs(maximum(abs.(Δi_diag)) - maximum(abs.(Δi_diag_prev)))
        push!(conv, (ΔΔ / maximum(Δi_diag)))
        push!(max_Δ, maximum(Δi_diag))
        if !isnothing(tol) && (ΔΔ / maximum(Δi_diag) <= tol) && (n > 1)
            # if this change is below a certain tolerance, stop iterating and return 
            return Δi_diag, conv
        end

        # if not, keep iterating 
        Δi_diag_prev = Δi_diag
    end

    return Δi_diag_prev, conv, max_Δ # the converged value for Δ
end

function compute_Δ(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real, periodic::Bool=true, niter::Int=100, tol::Union{Real,Nothing}=nothing, θ::Union{Real,Nothing})
    fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, periodic=periodic))

    # converge the BdG 
    Δi, conv, max_Δ = converge_BdG(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol, θ=θ, niter=niter, periodic=periodic, fsgap=fsgap)

    # the gap is the maximum value of Δi 
    #Δ = maximum(abs.(Δi))
end
