## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

# Parameters 
L = 11
t = 1
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
ϕx = 1.3
ϕy = 0.4
ϕz = 0.7
periodic = true
ndims = 3
T = 0

# We need to transform each of the eigenvectors into 2D space! 
function fourier_transform(m::ModelParams, u; minus=false)

    ndims, L = m.ndims, m.L

    # first, map back to 2D space 
    if ndims == 2
        N = L * L
        u = reshape(u, L, L)
    elseif ndims == 3
        N = L * L * L
        u = reshape(u, L, L, L)
    else
        println("$dims dimensions not implemented")
        return
    end

    # take FT along all spatial dimensions 
    if !minus # this is the expression for U_{q}
        uq = fft(u)
    else # this is the expression for U_{-q}
        uq = conj.(fft(conj.(u)))
    end

    # reshape it back and normalize. FFTW does not normalize!!
    uq = reshape(uq, N) ./ √N

    return uq
end

function IPR_momentum(m::ModelParams, u)
    # transform to momentum space
    uq = fourier_transform(m, u)
    u2 = uq .* conj.(uq)
    u4 = u2 .^ 2
    return real(sum(u4))
end

function IPR_real(u)
    u2 = u .* conj.(u)
    u4 = u2 .^ 2
    return real(sum(u4))
end

Js = LinRange(0, 3, 10)

dirpath = joinpath(@__DIR__, "hamiltonians")
mkpath(dirpath)

swave_IPR = zeros(2, length(Js))
dwave_IPR = zeros(2, length(Js))

for (Jᵢ, J) in enumerate(Js)
    m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)
    loadpath = joinpath(dirpath, "diagonal_$(J)J_$(L)L.h5")

    # calculate IPR in both real and momentum space for the susceptibility
    χ_swave, evec_swave = @time susceptibility_eigenvalue(m, T=T, symmetry="s-wave", calculate_dχdlogT=false, checkpointpath=loadpath, return_evec=true)
    # normalize the evec!! 
    evec_swave ./= norm(evec_swave)

    # χ_dwave, evec_dwave = @time susceptibility_eigenvalue(m, T=T, symmetry="d-wave", calculate_dχdlogT=false, checkpointpath=loadpath, return_evec=true)
    # Δ_dwave = to_7N_LGE_Δ(evec_dwave; L=m.L)
    # Δ_dwave = spatial_profile(m, Δ=Δ_dwave)
    # Δ_dwave = Δ_dwave[1, :, :, :]
    # reshape(Δ_dwave, L * L * L)

    swave_IPR[1, Jᵢ] = IPR_momentum(m, evec_swave)
    swave_IPR[2, Jᵢ] = IPR_real(evec_swave)
    # dwave_IPR[1, Jᵢ] = IPR_momentum(m, evec_dwave)
    # dwave_IPR[2, Jᵢ] = IPR_real(evec_dwave)

    # @assert 1 == 0
end

plot(Js, swave_IPR[1, :], c="blue", label=nothing)
plot!(Js, swave_IPR[2, :], c="red", label=nothing)
scatter!(Js, swave_IPR[1, :], label="Momentum IPR", c="blue")
scatter!(Js, swave_IPR[2, :], label="Real IPR", c="red")
title!("IPR of the maximum eigenvalue for $L×$L×$L system")
xlabel!("J")
ylabel!("IPR")