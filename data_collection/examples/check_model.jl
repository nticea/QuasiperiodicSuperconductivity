## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using CSV, DataFrames, Dates, Plots
include("../../src/BdG.jl")

# Parameters 
t = 1
J = 0
L = 7
Q = (√5 - 1) / 2
μ = 1e-8
θ = π / 7
V0 = 0
V1 = 0
ϕy = 2π * rand()
ϕz = 2π * rand()
ϕx = 2π * rand()
periodic = true
disorder = false
ndims = 3

m = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=ϕx, ϕy=ϕy, ϕz=ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
E, U = diagonalize_hamiltonian(m)

m̃ = ModelParams(L=L, t=t, Q=Q, μ=μ, θ=θ, ϕx=-ϕx, ϕy=-ϕy, ϕz=-ϕz, V0=V0, V1=V1, J=J, periodic=periodic, ndims=ndims, disorder=disorder)
Ẽ, Ũ = diagonalize_hamiltonian(m̃)

@assert Ẽ == E
@assert conj.(Ũ) == U

println("tests passed!")

# N = numsites(m)

# # We need to transform each of the eigenvectors into 2D space! 
# function fourier_transform_U(u; minus=false)
#     # first, map back to 2D space 
#     if ndims == 2
#         u = reshape(u, L, L)
#     elseif ndims == 3
#         u = reshape(u, L, L, L)
#     else
#         println("$dims dimensions not implemented")
#         return
#     end

#     # take FT along all spatial dimensions 
#     if !minus # this is the expression for U_{q}
#         uq = fft(u)
#     else # this is the expression for U_{-q}
#         uq = conj.(fft(conj.(u)))
#     end

#     # reshape it back and normalize. FFTW does not normalize!!
#     uq = reshape(uq, N) ./ √N

#     return uq
# end

# # perform a Fourier transform along the real space dim
# Uq = fourier_transform_U.(eachcol(U))
# Uq = hcat(Uq...)