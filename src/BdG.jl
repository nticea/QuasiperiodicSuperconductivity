using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations
using Distributions

include("../src/meanfield.jl")
include("../src/model.jl")

function finite_size_gap(; L::Int, t::Real, Q::Real, μ::Real, periodic::Bool=true, θ::Union{Real,Nothing}=nothing, ϕx::Real=0, ϕy::Real=0)
    H0 = noninteracting_hamiltonian(L=L, t=t, J=0, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
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

    @einsimd Δnew[i] := -V0 / 2 * U[i, n] * conj(V[i, n]) * tanh(E[n] / (2 * T))

    return Δnew, U, V, E
end

function rms(a, b)
    @assert size(a) == size(b)
    norm(a - b) / prod(size(a))
end

function converge_BdG(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool,
    V0::Real, θ::Union{Real,Nothing}, ϕx::Real=0, ϕy::Real=0, niter::Int=100,
    tol::Union{Real,Nothing}=nothing, fsgap::Union{Real,Nothing}=nothing, noise::Real=0)

    N = L * L

    # make the BdG equation matrix (fill in just block diagonals)
    M = zeros(2 * N, 2 * N)
    hij = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
    M[1:N, 1:N] .= hij
    M[(N+1):end, (N+1):end] .= -conj.(hij)

    if isnothing(fsgap)
        fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic))
    end

    # initial gues for Δ_i
    Δi = fsgap * Matrix(I, (N, N))
    Δi_diag_prev = diag(Δi)
    U_prev = zeros(N, 2 * N)
    V_prev = zeros(N, 2 * N)
    E_prev = zeros(2 * N)

    # iterate 
    conv = []
    max_Δ = []
    for n in 1:niter
        print(n, "-")
        Δi_diag, U, V, E = BdG_iteration(M, Δi; V0=V0, T=T) # perform one BdG iteration and get the new Δi 

        # add noise if desired (noise is set to 0 by default)
        d = Normal(0.0, maximum(Δi_diag) * noise)
        ε = rand.(d, size(Δi_diag)...)
        Δi_diag .+= ε

        Δi = diagm(Δi_diag) # make a matrix with Δi along the diagonals

        # calculate convergence information  
        push!(max_Δ, maximum(Δi_diag))

        Δi_diag_prev = Δi_diag
        U_prev = U
        V_prev = V
        E_prev = E
    end

    return Δi_diag_prev, U_prev, V_prev, E_prev, max_Δ # the converged value for Δ
end

function compute_Δ(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, periodic::Bool=true, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0)

    fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic))

    # converge the BdG 
    Δi, U, V, E, max_Δ = converge_BdG(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, fsgap=fsgap, noise=noise)

    return Δi, max_Δ
end

function BdG_coefficients(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, periodic::Bool=true, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0)

    fsgap = maximum(finite_size_gap(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic))

    # converge the BdG 
    Δi, U, V, E, max_Δ = converge_BdG(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, fsgap=fsgap, noise=noise)

    return U, V, E
end
