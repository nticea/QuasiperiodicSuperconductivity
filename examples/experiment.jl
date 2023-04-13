## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Plots
using LinearAlgebra

function V(kn, km)
    V(kn + km)
end

function V(k)
    cos(k)
    # 1 / 2 * (-2 + cos(2 * kx) + 2 * cos(ky) + cos(kx) * (2 + 4 * cos(ky)) + cos(2 * ky)) * sec(ky / 2)^2
end

niter = 100
N = 50
ks = LinRange(-2 * π, 2 * π, N)
err = []

global g0 = rand(N)
for _ in 1:niter
    for kn in ks
        for km in ks
            g = V(kn, km) ./ g0
            global g0 = copy(g)
        end
    end
end

plot(err)
