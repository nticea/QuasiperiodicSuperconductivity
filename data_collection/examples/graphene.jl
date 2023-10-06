## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Distributions, MCIntegration, Plots, Interpolations

function fermi(E::Real, T::Real; μ::Real)
    ξ = E - μ
    return 1 / (exp(ξ / T) + 1)
end

function make_interpolated_function(x, y, z)
    itp = LinearInterpolation((x, y), z)
    return itp
end

function plot_dos(B, dos_fcn; Bmax, num_LL, Δμ_offset)

    offset_fcn = zero_pt(Bmax=Bmax, num_LL=num_LL, numpts=numpts, dos_fcn=dos_fcn, Δμ_offset=Δμ_offset)
    offset = offset_fcn[B]

    Es = LinRange(-Elim, Elim, 500)
    d = []
    for E in Es
        push!(d, dos_fcn(E, B; num_LL=num_LL, Δμ_offset=Δμ_offset) - offset / 2)
    end

    return plot(Es, d, xlabel="E", ylabel="DOS", label=nothing)
end

function dos_monolayer(E, B; num_LL::Int=5, σ::Real=0.075, Δμ_offset::Real=0)
    E = E - Δμ_offset
    height = 4 * B
    dists = Normal{Float64}[]
    for L in 1:num_LL
        d1 = Normal(sqrt(L * B), σ)
        d2 = Normal(-sqrt(L * B), σ)
        push!(dists, d1, d2)
    end
    sum_dist = MixtureModel(dists)
    return pdf(sum_dist, E) * height * num_LL / 50
end

function dos_bernal(E, B; num_LL::Int=5, σ::Real=0.075, Δμ_offset::Real=0)
    E = E - Δμ_offset
    if E == 0
        height = 8 * B
    else
        height = 4 * B
    end
    dists = Normal{Float64}[]
    for L in 1:num_LL
        d1 = Normal(sqrt(L * (L - 1)) * B, σ)
        d2 = Normal(-sqrt(L * (L - 1)) * B, σ)
        push!(dists, d1, d2)
    end
    sum_dist = MixtureModel(dists)
    return pdf(sum_dist, E) * height * num_LL / 50
end

function dos_sum(E, B; num_LL::Int=5, σ::Real=0.2, Δμ_offset::Real=0)
    # we add the offset only to the bernal layer 
    db = dos_bernal(E, B, num_LL=num_LL, Δμ_offset=Δμ_offset)
    dm = dos_monolayer(E, B, num_LL=num_LL)

    return db + dm
end

function density_monolayer(; T::Real, Bmax::Real, num_LL::Int, numpts::Int, Δμ_offset::Real=0)
    return density(T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, dos_fcn=dos_monolayer, Δμ_offset=Δμ_offset)
end

function density_bernal(; T::Real, Bmax::Real, num_LL::Int, numpts::Int, Δμ_offset::Real=0)
    return density(T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, dos_fcn=dos_bernal, Δμ_offset=Δμ_offset)
end

function density_sum(; T::Real, Bmax::Real, num_LL::Int, numpts::Int, Δμ_offset::Real=0)
    return density(T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, dos_fcn=dos_sum, Δμ_offset=Δμ_offset)
end

function zero_pt(; Bmax::Real, num_LL::Int, numpts::Int, dos_fcn, Δμ_offset::Real=0)
    Bs = LinRange(0, Bmax, numpts)

    # integration variable 
    Evar = MCIntegration.Continuous(-Elim, Elim)

    # dofs 
    dof = [[1,] for _ in 1:length(Bs)]

    res = MCIntegration.integrate(var=Evar, dof=dof, print=-1, userdata=Bs) do intvar, c
        Bs = c.userdata

        function density_calc(B)
            E = intvar[1]
            return fermi(E, T, μ=0) * dos_fcn(E, B, num_LL=num_LL, Δμ_offset=Δμ_offset)
        end

        return [density_calc(B) for B in Bs]
    end

    ns = res.mean
    itp = LinearInterpolation(Bs, ns / Elim)

    return itp
end

function density(; T::Real, Bmax::Real, num_LL::Int, numpts::Int, dos_fcn, Δμ_offset::Real=0)
    Bs = LinRange(0, Bmax, numpts)
    μs = LinRange(-Elim, Elim, 2 * numpts)

    Bμs = [[B, μ] for B in Bs, μ in μs]
    Bμs = reshape(Bμs, prod(size(Bμs)))

    # integration variable 
    Evar = MCIntegration.Continuous(-Elim, Elim)

    # dofs 
    dof = [[1,] for _ in 1:length(Bμs)]

    # subtract away the mean value 
    offset_fcn = zero_pt(Bmax=Bmax, num_LL=num_LL, numpts=numpts, dos_fcn=dos_fcn, Δμ_offset=Δμ_offset)

    res = MCIntegration.integrate(var=Evar, dof=dof, print=-1, userdata=Bμs) do intvar, c
        Bμs = c.userdata

        function density_calc(Bμ)
            B, μ = Bμ
            offset = offset_fcn[B]
            E = intvar[1]
            return fermi(E, T, μ=μ) * dos_fcn(E, B, num_LL=num_LL, Δμ_offset=Δμ_offset) - offset / 2
        end

        return [density_calc(Bμ) for Bμ in Bμs]
    end

    n = res.mean
    n = reshape(res.mean, length(Bs), length(μs))
    return n
end

function dos_nB(kind::String; T::Real, Bmax::Real, num_LL::Int, numpts::Int, Δμ_offset::Real=0)
    if kind == "bernal"
        dos_fcn = dos_bernal
        n = density_bernal(T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, Δμ_offset=Δμ_offset)
    elseif kind == "monolayer"
        dos_fcn = dos_monolayer
        n = density_monolayer(T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, Δμ_offset=Δμ_offset)
    elseif kind == "sum"
        dos_fcn = dos_sum
        n = density_sum(T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, Δμ_offset=Δμ_offset)
    else
        @error "$kind not recognized"
        return
    end

    Bs = LinRange(0, Bmax, numpts)
    μs = LinRange(-Elim, Elim, 2 * numpts)

    # subtract away the offset 
    # offset_fcn = zero_pt(Bmax=Bmax, num_LL=num_LL, numpts=numpts, dos_fcn=dos_fcn, Δμ_offset=Δμ_offset)
    # offset = offset_fcn[Bs]

    # density of states result for all Bs 
    dos = [dos_fcn(μ, B, num_LL=num_LL, Δμ_offset=Δμ_offset) for B in Bs, μ in μs]

    ns_result = range(-1, 1, length=500)
    dos_result = zeros(length(Bs), length(ns_result))
    for (b, B) in enumerate(Bs)
        print("$b-")
        nB = n[b, :]
        dosB = dos[b, :]

        sortidx = sortperm(nB)
        nB = nB[sortidx]
        dosB = dosB[sortidx]

        Interpolations.deduplicate_knots!(nB, move_knots=true)

        dos_itp = LinearInterpolation(nB, dosB)
        for (i, nᵢ) in enumerate(ns_result)
            try
                dos_result[b, i] = dos_itp[nᵢ]
            catch
                dos_result[b, i] = NaN
            end
        end

    end
    return dos_result
end

# external variables
Elim = 5.4
T = 1e-3
Bmax = 5
num_LL = 10
numpts = 100
Δμ_offset = 1

plotsum = plot_dos(3, dos_sum; Bmax, num_LL, Δμ_offset)
plotb = plot_dos(3, dos_bernal; Bmax, num_LL, Δμ_offset)
plotm = plot_dos(3, dos_monolayer; Bmax, num_LL, Δμ_offset=0)

dos_nB_bernal = dos_nB("bernal", T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, Δμ_offset=Δμ_offset)
dos_nB_monolayer = dos_nB("monolayer", T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, Δμ_offset=0)
dos_nB_sum = dos_nB("sum", T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, Δμ_offset=Δμ_offset)

Bs = LinRange(0, Bmax, numpts)
ns = range(-1, 1, length=500)
hb = heatmap(ns, Bs, dos_nB_bernal, xlabel="n", ylabel="B", title="Bernal")
hm = heatmap(ns, Bs, dos_nB_monolayer, xlabel="n", ylabel="B", title="Monolayer")
hs = heatmap(ns, Bs, dos_nB_sum, xlabel="n", ylabel="B", title="Bernal + Monolayer", c=reverse(cgrad(:magma)))

# n_fcn = density_monolayer(T=T, Bmax=Bmax, num_LL=num_LL, numpts=numpts, Δμ_offset=Δμ_offset)

# B = 3
# μs = range(-Elim, Elim, length=500)
# n = [n_fcn(B, μ) for μ in μs]
# plot(μs, n, xlabel="μ", ylabel="n", legend=nothing, title="n(μ) at B=$B")

# x = dos_nB_bernal[:, 1] # B
# y = dos_nB_bernal[:, 2] # n
# z = dos_nB_bernal[:, 3] # DOS

# # Define the boundaries and size of grid cells
# Δx = 0.01  # Adjust as needed
# Δy = 0.01  # Adjust as needed

# # Create arrays for bin edges
# x_bins = minimum(x):Δx:maximum(x)
# y_bins = minimum(y):Δy:maximum(y)

# # Initialize an array to store the binned data
# binned_data = zeros(Float64, length(x_bins) - 1, length(y_bins) - 1)

# # Loop through data points and bin the values
# for i = 1:length(x)
#     x_val = x[i]
#     y_val = y[i]
#     z_val = z[i]

#     # Find the bin indices for x and y
#     x_bin = searchsortedfirst(x_bins, x_val)
#     y_bin = searchsortedfirst(y_bins, y_val)

#     # Add the z value to the corresponding bin
#     if 1 ≤ x_bin ≤ length(x_bins) - 1 && 1 ≤ y_bin ≤ length(y_bins) - 1
#         binned_data[x_bin, y_bin] += z_val
#     end
# end

# heatmap(binned_data)

## THIS WORKS HERE ## 

# dos_fcn = dos_bernal
# Bs = LinRange(0, Bmax, 50)
# μs = LinRange(-Elim, Elim, 50)

# Bμs = [[B, μ] for B in Bs, μ in μs]
# Bμs = reshape(Bμs, prod(size(Bμs)))

# # integration variable 
# Evar = MCIntegration.Continuous(-Elim, Elim)

# # dofs 
# dof = [[1,] for _ in 1:length(Bμs)]

# # subtract away the mean value 
# offset_fcn = zero_pt(Bmax=Bmax, num_LL=num_LL, numpts=numpts, dos_fcn=dos_fcn, Δμ_offset=Δμ_offset)

# res = MCIntegration.integrate(var=Evar, dof=dof, print=-1, userdata=Bμs) do intvar, c
#     Bμs = c.userdata

#     function density_calc(Bμ)
#         B, μ = Bμ
#         offset = offset_fcn[B]
#         E = intvar[1]
#         return fermi(E, T, μ=μ) * dos_fcn(E, B, num_LL=num_LL, Δμ_offset=Δμ_offset) - offset / 2
#     end

#     return [density_calc(Bμ) for Bμ in Bμs]
# end

# n = res.mean
# n = reshape(res.mean, length(Bs), length(μs))

# offset_fcn = zero_pt(Bmax=Bmax, num_LL=num_LL, numpts=numpts, dos_fcn=dos_fcn, Δμ_offset=Δμ_offset)
# offset = offset_fcn[Bs]
# dos = [dos_fcn(μ, B, num_LL=num_LL, Δμ_offset=Δμ_offset) for B in Bs, μ in μs]

# Bidx = Bmax
# plot(dos[Bidx, :])
# plot!(n[Bidx, :])

# plot(n[Bidx, :], dos[Bidx, :])