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

function BdG_iteration_swave(M::Matrix{Float64}, Δi; V0::Real, T::Real)
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

function initialize_Δ_swave(m::ModelParams; Δ_init=nothing)
    N = numsites(m)
    if isnothing(Δ_init)
        Δi = Matrix(I, (N, N))
        return Δi .* 0.05
    else
        return ΔLGE_to_ΔBdG(m, Δ=Δ_init)
    end
end

function converge_BdG_swave(m::ModelParams; T::Real, Δ_init=nothing, niter::Int=500, tol::Real=1e-12, scratchpath=nothing)
    N = numsites(m)

    # make the BdG equation matrix (fill in just block diagonals)
    M = zeros(2 * N, 2 * N)
    hij = noninteracting_hamiltonian(m)
    M[1:N, 1:N] .= hij
    M[(N+1):end, (N+1):end] .= -conj.(hij)

    # initial gues for Δ_i
    if !isnothing(scratchpath)
        df = load_BdG_checkpoint(m, scratchpath, symmetry)
        if size(df)[1] > 0
            nstart = df.n[1]
            Δ_init = df.Δ[1]
            println("Loading checkpointed Δ at iteration $nstart")
        else
            nstart = 1
        end
    end

    Δi = initialize_Δ_swave(m, Δ_init=Δ_init)
    Δi_prev = copy(Δi)
    U, V, E = nothing, nothing, nothing

    # iterate 
    max_Δ = []
    for n in nstart:niter
        print(n, "-")

        # perform one BdG iteration and get the new Δi 
        Δi_diag, U, V, E = BdG_iteration_swave(M, Δi; V0=V0, T=T)

        # make a matrix with Δi along the diagonals
        Δi = diagm(Δi_diag)

        # check for convergence 
        ΔΔ = norm(Δi .- Δi_prev)
        if ΔΔ <= tol
            Δij = ΔBdG_to_ΔLGE_flat(m, Δ=Δi)
            return Δij, U, V, E, max_Δ
        end

        # track convergence history  
        push!(max_Δ, maximum(Δi))

        Δi_prev = Δi

        if !isnothing(scratchpath)
            # checkpoint 
            Δij = ΔBdG_to_ΔLGE_flat(m, Δ=Δij)
            stamp = "BdG_$(m.ndims)D_$(m.L)L_$(m.J)J_$(round(m.θ, digits=3))theta_$(round(m.Q,digits=3))Q_$(m.V0)V0_$(m.V1)V1.csv"
            BdG_checkpoint = DataFrame(L=[m.L], t=[m.t], μ=[m.μ], J=[m.J],
                Q=[m.Q], θ=[m.θ], ϕx=[m.ϕx], ϕy=[m.ϕy], ϕz=[m.ϕz], V0=[m.V0], V1=[m.V1],
                ndims=[m.ndims], periodic=[m.periodic],
                disorder=[m.disorder], Δ=[Δi], n=[n])
            fname = joinpath(scratchpath, stamp)
            if isfile(fname)
                CSV.write(fname, BdG_checkpoint, append=true)
            else
                CSV.write(fname, BdG_checkpoint, append=false)
            end
            flush(stdout)
        end
    end

    Δij = ΔBdG_to_ΔLGE_flat(m, Δ=Δi)
    return Δij, U, V, E, max_Δ
end

function converge_BdG_dwave(m::ModelParams; T::Real, Δ_init=nothing, niter::Int=500, tol::Real=1e-12, scratchpath=nothing)
    converge_BdG_nonswave(m; T=T, Δ_init=Δ_init, niter=niter, tol=tol, scratchpath=scratchpath, symmetry="d-wave")
end

function converge_BdG_pwave(m::ModelParams; T::Real, Δ_init=nothing, niter::Int=500, tol::Real=1e-12, scratchpath=nothing)
    converge_BdG_nonswave(m; T=T, Δ_init=Δ_init, niter=niter, tol=tol, scratchpath=scratchpath, symmetry="p-wave")
end

function converge_BdG_nonswave(m::ModelParams; T::Real, Δ_init=nothing, niter::Int=500, tol::Real=1e-12, scratchpath=nothing, symmetry::String)

    # size of the BdG matrix is 2N × 2N
    N = numsites(m)

    # make the BdG matrix (fill in just the diagonals with H0 for now)
    M = zeros(2 * N, 2 * N)
    hij = noninteracting_hamiltonian(m)
    M[1:N, 1:N] .= hij
    M[(N+1):end, (N+1):end] .= -conj.(hij)

    # make the interaction matrix 
    Vij = make_interaction(m)

    # initial guess for Δ_i
    if !isnothing(scratchpath)
        df = load_BdG_checkpoint(m, scratchpath, symmetry)
        if size(df)[1] > 0
            nstart = df.n[1]
            Δ_init = df.Δ[1]
            println("Loading checkpointed Δ at iteration $nstart")
        else
            nstart = 1
        end
    else
        nstart = 1
    end

    # make an initial guess for Δ
    if symmetry == "d-wave"
        Δij_prev1 = initialize_Δ_dwave(m, Δ_init=Δ_init)
    elseif symmetry == "p-wave"
        Δij_prev1 = initialize_Δ_pwave(m, Δ_init=Δ_init)
    else
        println("only p- and d-wave, sorry")
        return
    end
    Δij_prev2 = copy(Δij_prev1)

    U, V, E = nothing, nothing, nothing

    ## ITERATE ## 
    hist = [] # keep track of convergence history
    for n in nstart:niter
        print(n, "-")

        # take the average of the previous 2 iterations 
        Δij_prev = (Δij_prev1 + Δij_prev2) ./ 2

        # Compute the BdG iteration 
        if symmetry == "d-wave"
            Δij, U, V, E = BdG_iteration_dwave(M, Δij_prev; Vij=Vij, T=T) # perform one BdG iteration and get the new Δij 
        elseif symmetry == "p-wave"
            Δij, U, V, E = BdG_iteration_pwave(M, Δij_prev; Vij=Vij, T=T) # perform one BdG iteration and get the new Δij 
        else
            println("only p- and d-wave, sorry")
            return
        end

        # check for convergence 
        ΔΔ = norm(Δij .- Δij_prev)
        if ΔΔ <= tol
            Δij = ΔBdG_to_ΔLGE_flat(m, Δ=Δij)
            return Δij, U, V, E, hist
        end

        # keep track of convergence over time 
        push!(hist, ΔΔ)

        # Update the history 
        Δij_prev2 = copy(Δij_prev1)
        Δij_prev1 = copy(Δij)

        if !isnothing(scratchpath)
            # checkpoint 
            Δij = ΔBdG_to_ΔLGE_flat(m, Δ=Δij)
            stamp = "BdG_$(m.ndims)D_$(m.L)L_$(m.J)J_$(round(m.θ, digits=3))theta_$(round(m.Q,digits=3))Q_$(m.V0)V0_$(m.V1)V1.csv"
            BdG_checkpoint = DataFrame(L=[m.L], t=[m.t], μ=[m.μ], J=[m.J],
                Q=[m.Q], θ=[m.θ], ϕx=[m.ϕx], ϕy=[m.ϕy], ϕz=[m.ϕz], V0=[m.V0], V1=[m.V1],
                ndims=[m.ndims], periodic=[m.periodic],
                disorder=[m.disorder], Δ=[Δij], n=[n], symmetry=[symmetry])
            fname = joinpath(scratchpath, stamp)
            if isfile(fname)
                CSV.write(fname, BdG_checkpoint, append=true)
            else
                CSV.write(fname, BdG_checkpoint, append=false)
            end
            flush(stdout)
        end
    end

    Δij = (Δij_prev1 + Δij_prev2) ./ 2
    Δij = ΔBdG_to_ΔLGE_flat(m, Δ=Δij)

    return Δij, U, V, E, hist
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


function BdG_iteration_pwave(M::Matrix{Float64}, Δij; Vij, T::Real)
    M = copy(M)
    N = size(Δij)[1]

    # Fill in the off-diagonal blocks with Δ
    M[1:N, (N+1):end] .= Δij
    M[(N+1):end, 1:N] .= -conj.(Δij) # p-wave has neg sign out front

    # diagonalize this matrix to get the eigenvalues 
    E, UV = eigen(Hermitian(M))

    # Extract the u and the v coefficients
    U = UV[1:N, :]
    V = UV[(N+1):end, :]

    # compute the updated Δij
    # @einsimd Δnew[i, j] := Vij[i, j] / 4 * (U[i, n] * conj(V[j, n]) - U[j, n] * conj(V[i, n])) * tanh(E[n] / (2 * T))
    fs = fermi.(E, T)
    fsneg = fermi.(-E, T)
    @einsimd Δnew[i, j] := Vij[i, j] * (U[i, n] * conj(V[j, n]) * fsneg[n] + U[j, n] * conj(V[i, n]) * fs[n])

    return Δnew, U, V, E
end

function initialize_Δ_dwave(m::ModelParams; Δ_init=nothing)
    N = numsites(m)
    if isnothing(Δ_init)
        d = Uniform(0, 0.05)
        return rand.(d, N, N)
    else
        return ΔLGE_to_ΔBdG(m, Δ=Δ_init)
    end
end

function initialize_Δ_pwave(m::ModelParams; Δ_init=nothing)
    N = numsites(m)
    if isnothing(Δ_init)
        d = Uniform(0, 0.05)
        mat = rand.(d, N, N)
        mat[diagind(mat)] .= 0
        return mat
    else
        return ΔLGE_to_ΔBdG(m, Δ=Δ_init)
    end
end

function make_interaction(m::ModelParams)
    L, V0, V1, periodic = m.L, m.V0, m.V1, m.periodic

    # construct the adjacency matrix for the nearest-neighbours 
    if ndims == 2
        g = Graphs.SimpleGraphs.grid((L, L), periodic=periodic)
    elseif ndims == 3
        g = Graphs.SimpleGraphs.grid((L, L, L), periodic=periodic)
    else
        println("ndims=$ndims not supported (yet?)")
        return
    end
    Vij = Graphs.LinAlg.adjacency_matrix(g)
    Vij = Matrix(Vij) # temporarily make it a dense matrix 
    Vij = convert(Matrix{Float64}, Vij)

    # the off-diagonals are multiplied by V1 
    Vij .*= V1

    # the diagonals are multiplied by V0 
    Vij[diagind(Vij)] .= V0

    return Vij
end

function compute_Δ(m; T::Real, Δ_init=nothing, niter::Int=100, tol::Union{Real,Nothing}=nothing, symmetry::String)

    # converge the BdG 
    if m.V1 == 0
        Δi, _, _, _, max_Δ = converge_BdG_swave(m, T=T, niter=niter, tol=tol, Δ_init=Δ_init)
    else
        if symmetry == "d-wave"
            Δi, _, _, _, max_Δ = converge_BdG_dwave(m, T=T, niter=niter, tol=tol, Δ_init=Δ_init)
        elseif symmetry == "p-wave"
            Δi, _, _, _, max_Δ = converge_BdG_pwave(m, T=T, niter=niter, tol=tol, Δ_init=Δ_init)
        else
            println("sucks.")
            return
        end
    end

    # normalize 
    Δi ./= norm(Δi)

    return Δi, max_Δ
end

function BdG_coefficients(m; T::Real, Δ_init=nothing, niter::Int=100, tol::Union{Real,Nothing}=nothing, scratchpath=nothing, symmetry::String="d-wave")

    # converge the BdG 
    if m.V1 == 0
        Δi, U, V, E, _ = converge_BdG_swave(m, T=T, niter=niter, tol=tol, Δ_init=Δ_init, scratchpath=scratchpath)
    else
        if symmetry == "d-wave"
            Δi, U, V, E, _ = converge_BdG_dwave(m, T=T, niter=niter, tol=tol, Δ_init=Δ_init, scratchpath=scratchpath)
        elseif symmetry == "p-wave"
            Δi, U, V, E, _ = converge_BdG_pwave(m, T=T, niter=niter, tol=tol, Δ_init=Δ_init, scratchpath=scratchpath)
        else
            println("sucks.")
            return
        end
    end

    # normalize 
    Δi ./= norm(Δi)

    return U, V, E, Δi
end

function load_BdG_checkpoint(m::ModelParams, dirname, symmetry)
    # read files 
    files = readdir(joinpath(@__DIR__, dirname))
    df = DataFrame(L=[], Δ=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[],
        ϕz=[], ndims=[], V0=[], V1=[], t=[], μ=[], periodic=[],
        disorder=[], n=[], symmetry=[])
    for f in files
        if endswith(f, ".csv") && contains(f, "BdG")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, dirname, f)))
                dfi = convert_df_arrays(dfi, "Δ")
                append!(df, dfi)
            catch e
                @show e
                println("file $f not valid...?")
            end
        end
    end
    dropmissing!(df)
    df = df[(df.L.==m.L).&&(df.J.==m.J).&&(df.Q.==m.Q).&&(df.θ.==m.θ).&&(df.ϕx.==m.ϕx).&&(df.ϕy.==m.ϕy).&&(df.ϕz.==m.ϕz).&&(df.ndims.==m.ndims).&&(df.V0.==m.V0).&&(df.V1.==m.V1).&&(df.t.==m.t).&&(df.μ.==m.μ).&&(df.periodic.==m.periodic).&&(df.disorder.==m.disorder).&&(df.symmetry.==symmetry), :]
    sort!(df, [:n], rev=true)
    return df
end