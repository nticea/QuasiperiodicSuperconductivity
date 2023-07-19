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
include("../src/results.jl")

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

function initialize_Δ_swave(; L::Int, Δ_init=nothing)
    # N = L * L
    # if isnothing(Δ_init)
    #     Δi = Matrix(I, (N, N))
    #     d = Uniform(0, 0.05)
    #     Δi[diagind(Δi)] = rand.(d, L * L)
    #     return Δi
    # else
    #     @assert size(Δ_init) == (5 * L * L,)
    #     return ΔLGE_to_ΔBdG(Δ_init, L=L)
    # end
    N = L * L
    Δi = Matrix(I, (N, N))
    # d = Uniform(0, 0.05)
    # Δi[diagind(Δi)] = rand.(d, L * L)
    return Δi .* 0.05
end

function converge_BdG(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool,
    V0::Real, θ::Union{Real,Nothing}, ϕx::Real=0, ϕy::Real=0, niter::Int=100,
    tol::Union{Real,Nothing}=nothing, fsgap::Union{Real,Nothing}=nothing, noise::Real=0,
    Δ_init=nothing)

    N = L * L

    # make the BdG equation matrix (fill in just block diagonals)
    M = zeros(2 * N, 2 * N)
    hij = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, periodic=periodic)
    M[1:N, 1:N] .= hij
    M[(N+1):end, (N+1):end] .= -conj.(hij)

    # initial gues for Δ_i
    Δi = initialize_Δ_swave(L=L, Δ_init=Δ_init)

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

    # converge the BdG 
    Δi, U, V, E, max_Δ = converge_BdG(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise)

    return Δi, max_Δ
end

function BdG_coefficients_swave(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, periodic::Bool=true, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0, Δ_init=nothing)

    # converge the BdG 
    Δi, U, V, E, _ = converge_BdG(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol,
        θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise, Δ_init=Δ_init)

    return U, V, E, Δi
end
