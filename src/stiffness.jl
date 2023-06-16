include("BdG.jl")
include("BdG_dwave.jl")

function superfluid_stiffness(T; L::Int, t::Real, J::Real, Q::Real, μ::Real, periodic::Bool, V0::Real, V1::Real, θ::Union{Real,Nothing},
    ϕx::Real=0, ϕy::Real=0, niter::Int=100, tol::Union{Real,Nothing}=nothing, noise::Real=0)

    if V1 == 0
        U, V, E = BdG_coefficients(T, L=L, t=t, J=J, Q=Q, μ=μ, V0=V0, tol=tol, θ=θ, ϕx=ϕx, ϕy=ϕy, niter=niter, periodic=periodic, noise=noise)
    else
        @error "V1 ≂̸ 0 Not implemented yet"
    end

    f = fermi.(E; T=T)

    N, _ = size(U)
    Ds = []
    for i in 1:N
        nn = [nearest_neighbours(i, L=L)...] # get the nearest neighbours
        D_i = []
        for j in nn # order of nn: right, up, left, down 
            Kx = kinetic_term(i, j; U=U, V=V, f=f, t=t)
            Πxx = current_current_term(i, j; L=L, U=U, V=V, f=f, t=t)

            push!(D_i, -Kx + Πxx)
        end

        push!(Ds, D_i)
    end

end

function kinetic_term(i, j; U, V, f, t)
    @error "I might be missing a factor of 2 from the spin σ"
    N, _ = size(U)

    U_i = U[i, :]
    U_j = U[j, :]
    Uconj_i = conj.(U)[i, :]
    Uconj_j = conj.(U)[j, :]

    V_i = V[i, :]
    V_j = V[j, :]
    Vconj_i = conj.(V)[i, :]
    Vconj_j = conj.(V)[j, :]

    @einsimd term1 := f[n] * (Uconj_j[n] * U_i[n] + Uconj_i[n] * U_j[n])
    @einsimd term2 := (1 - f[n]) * (Vconj_j[n] * V_i[n] + Vconj_i[n] * V_j[n])

    return -t / N * (term1 + term2)
end

function current_current_term(i, j; L, U, V, f, t)
    x, y = site_to_coordinate(i, L=L)

    exp(-1im * ())
end

function Aq(q)

end

function Dq(q)

end

function fermi(ε; T)
    1 / (1 + exp(ε / T))
end