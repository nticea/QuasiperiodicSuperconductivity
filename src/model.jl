using LinearAlgebra
using SparseArrays
using Graphs
using ArnoldiMethod
using TriangularIndices
using Tullio
using Einsum
using Interpolations

include("../src/results.jl")

function expspace(start, stop, length)
    exp10.(range(start, stop, length=length))
end

function λmax(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, V0::Real)
    # Construct the non-interacting Hamiltonian matrix
    H0 = noninteracting_hamiltonian(L=L, t=t, J=J, Q=Q, μ=μ)
    #heatmap(Matrix(H0),yflip=true,clims=(-maximum(abs.(H0)),maximum(abs.(H0))),cmap=:bwr)

    # Diagonalize this Hamiltonian
    E, U = diagonalize_hamiltonian(H0)

    # Construct the pairfield susceptibility
    χ = pairfield_singlet(T, E=E, U=U)

    # Construct M (for s-wave, all we do is multiply χ by +V0)
    M = make_M(χ, V0)

    # Calculate Tc by finding the eigenvalues of M
    λs = diagonalize_M(M)

    return λs[1]
end

function noninteracting_hamiltonian(; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool=true)
    # construct the kinetic part 
    Ht = square_lattice_kinetic(L=L, t=t, periodic=periodic)

    # interaction
    Hint = zeros(L * L) # this is just the diagonal 
    for x in 1:L
        for y in 1:L
            n = coordinate_to_site(x, y, L=L)
            U_xy = aubry_andre(x, y, J=J, Q=Q)
            Hint[n] = -(U_xy + μ)
        end
    end
    Hint = spdiagm(Hint)

    return Ht + Hint
end

function square_lattice_kinetic(; L::Int, t::Real, periodic::Bool=true)
    g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    H = Graphs.LinAlg.adjacency_matrix(g)
    return -t .* H
end

function coordinate_to_site(x::Int, y::Int; L::Int)
    (x - 1) * L + y
end

function site_to_coordinate(r; L::Int)
    x = floor(Int, r / L) + 1
    y = r % L
    return x, y
end

function aubry_andre(x, y; J::Real, Q::Real)
    J * (cos(2 * π * Q * (x + y)) - cos(2 * π * Q * (x - y)))
end

function fermi(ε::Real, T::Real)
    # note: the chemical potential was subtracted away in the Hamiltonian 
    1 / (exp(ε / T) + 1)
end

function diagonalize_hamiltonian(H)
    # we must compute all eigenvalues
    vals, vecs = eigen(Hermitian(Matrix(H)))
    return vals, vecs
end

function diagonalize_M(M)
    decomp, _ = partialschur(M, nev=1, tol=1e-6, which=LM())
    return decomp.R
end

function pairfield_singlet(T::Real; E, U)
    N = size(U)[1]
    χ = zeros(N, N)

    Uconj = conj.(U)

    # make the prefactor
    fs = fermi.(E, T)
    fnm = zeros(N, N)
    Enm = zeros(N, N)
    for i in 1:N
        fnm[:, i] = fs .+ fs[i]
        Enm[:, i] = E .+ E[i]
    end
    P = (1 .- fnm) ./ Enm

    @einsimd PUU[r, n, m] := Uconj[r, n] * P[n, m] * Uconj[r, m]
    @einsimd UU[m, n, rprime] := U[rprime, n] * U[rprime, m]
    PUU = reshape(PUU, N, N * N)
    UU = reshape(UU, N * N, N)
    χ = PUU * UU

    return χ
end

function make_M(χ, V0)
    return χ * V0 # scale by V0
end

function find_Tc(results::Results)
    L, λs, Js, V0s, Ts = results.L, results.λs, results.Js, results.V0s, results.Ts

    # do interpolations 
    Tcs = zeros(length(Js), length(V0s))
    for k in 1:length(Js)
        for j in 1:length(V0s)
            knots = reverse(λs[k, j, :])
            Interpolations.deduplicate_knots!(knots, move_knots=true)
            try
                interp_linear = linear_interpolation(knots, reverse(Ts))
                Tcs[k, j] = interp_linear(1)
            catch e
                Tcs[k, j] = NaN
            end
        end
    end

    return Tcs
end

function plot_Tcs(results::Results)
    L, λs, Js, V0s, Ts = results.L, results.λs, results.Js, results.V0s, results.Ts
    Tcs = find_Tc(results)

    p2 = plot()
    cmap = cgrad(:Set1_9, length(V0s), categorical=true)
    for (k, J) in enumerate(Js)
        plot!(p2, V0s, Tcs[k, :], xaxis=:log10, yaxis=:log10, c=cmap[k], label=nothing)
        scatter!(p2, V0s, Tcs[k, :], xaxis=:log10, yaxis=:log10, c=cmap[k], label="J=$(J)")
    end

    title!(p2, "Transition temperature for $(L)x$(L) square lattice")
    xlabel!(p2, "V")
    ylabel!(p2, "Tc")
end