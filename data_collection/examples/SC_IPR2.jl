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
ϕx = 0.9
ϕy = 0.8
ϕz = 0
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

Js_full = LinRange(0, 3, 10)
Js = [Js_full[2], Js_full[9], 7]

dirpath = joinpath(@__DIR__, "hamiltonians")
mkpath(dirpath)

J = Js[2]
m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=-1, V1=-1, J=J, periodic=periodic, ndims=ndims)
loadpath = joinpath(dirpath, "diagonal_$(J)J_$(L)L.h5")

# calculate IPR in both real and momentum space for the susceptibility
χ_swave, evec_swave = @time susceptibility_eigenvalue(m, T=T, symmetry="s-wave", calculate_dχdlogT=false, checkpointpath=loadpath, return_evec=true)
# normalize the evec!! 
Δ = spatial_profile(m, Δ=evec_swave)
plot_spatial_profile(m, Δ=Δ)