using DataFrames
using HDF5
using Statistics

include("IndexMatrix.jl")
include("model.jl")

function update_results!(m::ModelParams, df::DataFrame; T, Δ, λ)
    L, t, Q, μ, J, θ, ϕx, ϕy, ϕz, V0, V1, ndims, periodic = m.L, m.t, m.Q, m.μ, m.J, m.θ, m.ϕx, m.ϕy, m.ϕz, m.V0, m.V1, m.ndims, m.periodic
    df2 = DataFrame(L=[L], t=[t], Q=[Q], μ=[μ], J=[J], θ=[θ], ϕx=[ϕx],
        ϕy=[ϕy], ϕz=[ϕz], V0=[V0], V1=[V1], ndims=[ndims], periodic=[periodic],
        T=[T], Δ=[Δ], λ=[λ])
    append!(df, df2)
end

function update_results!(m::ModelParams, df::DataFrame; λ, T, Δ, K=zeros(4), Π=zeros(4))
    L, t, Q, μ, J, θ, ϕx, ϕy, ϕz, V0, V1, ndims, periodic = m.L, m.t, m.Q, m.μ, m.J, m.θ, m.ϕx, m.ϕy, m.ϕz, m.V0, m.V1, m.ndims, m.periodic
    df2 = DataFrame(L=[L], t=[t], Q=[Q], μ=[μ], J=[J], θ=[θ], ϕx=[ϕx],
        ϕy=[ϕy], ϕz=[ϕz], V0=[V0], V1=[V1], ndims=[ndims], periodic=[periodic],
        T=[T], Δ=[Δ], λ=[λ], K=[K], Π=[Π])
    append!(df, df2)
end

function already_calculated(m::ModelParams, df::DataFrame; T)
    L, t, Q, μ, J, θ, ϕx, ϕy, ϕz, V0, V1, ndims, periodic = m.L, m.t, m.Q, m.μ, m.J, m.θ, m.ϕx, m.ϕy, m.ϕz, m.V0, m.V1, m.ndims, m.periodic
    sub = df[(df.L.==L).&(df.t.==t).&(df.μ.==μ).&(df.J.==J).&(df.Q.==Q).&(df.V0.==V0).&(df.V1.==V1).&(df.θ.==θ).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy).&(df.ϕz.==ϕz).&(df.T.==T).&(df.ndims.==ndims).&(df.periodic.==periodic), :]
    return size(sub)[1] > 0
end

# same as above, but without T (this is for finding Tc using LGE method)
function already_calculated(m::ModelParams, df::DataFrame)
    L, μ, J, θ, ϕx, ϕy, ϕz, V0, V1, ndims = m.L, m.μ, m.J, m.θ, m.ϕx, m.ϕy, m.ϕz, m.V0, m.V1, m.ndims
    sub = df[(df.L.==L).&(df.t.==t).&(df.μ.==μ).&(df.J.==J).&(df.Q.==Q).&(df.V0.==V0).&(df.V1.==V1).&(df.θ.==θ).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy).&(df.ϕz.==ϕz).&(df.ndims.==ndims).&(df.periodic.==periodic), :]
    return size(sub)[1] > 0
end

function convert_df_arrays(df::DataFrame, col_name::String, delims=r"[,; ]")
    all_arrs = []
    for (i, arr) in enumerate(df[!, Symbol(col_name)])
        cleaned_string = replace(arr, r"Any|\[|\]" => "")
        # Split the string by multiple delimiters
        result = split(cleaned_string, delims)
        new_array = result[result.!=""]
        arrs_new = parse.(Float64, new_array)
        push!(all_arrs, arrs_new)
    end
    dfcut = copy(df)#df[:, collect(1:size(df)[2]-1)]
    dfcut[!, Symbol(col_name)] = all_arrs

    return dfcut
end

function find_Tc_from_dataframe(df::DataFrame; interp_value::Real=1)
    println("TODO: Implement this using grouped dataframes")
    Tc_df = DataFrame([name => [] for name in names(df)])
    grouped = groupby(df, select(df, Not(:T)))
    for dfsub in grouped
        λs, Ts = dfsub.λ, dfsub.T
        # Compute the interpolated Tc for this (J,V0) pair
        idxlist = sortperm(λs)
        knots = λs[idxlist]
        if length(knots) > 1
            Interpolations.deduplicate_knots!(knots, move_knots=true)
            df2 = DataFrame([name => [dfsub[name][1]] for name in names(dfsub)])
            try
                interp_linear = linear_interpolation(knots, Ts[idxlist])
                df2.Tc = interp_linear(interp_value)
                append!(Tc_df, df2)
            catch e
                Tc = NaN
                df2.Tc = NaN
                append!(Tc_df, df2)
            end
        end
    end
    return Tc_df
end

function find_Tc_from_LGE_df(df)
    return find_Tc_from_dataframe(df, interp_value=1)
end

function find_Tc_from_BdG_df(df)
    return find_Tc_from_dataframe(df, interp_value=1e-5)
end

function θ_to_π(θ)
    if isnothing(θ)
        return "untilted"
    end
    for n in 1:10
        if isapprox(n * θ, π)
            return "π/$n"
        end
    end
    return θ
end

function colour_phase(x1::Int, x2::Int, x3::Int; all_evs, numpts::Int=10)
    if Base.ndims(all_evs) == 3
        val = all_evs[x1, x2, x3]
    elseif Base.ndims(all_evs) == 2
        val = all_evs[x2, x3]
    end
    if val < 0
        return "blue"
    else
        return "red"
    end
end

function spatial_profile(m::ModelParams; Δ)
    L, ndims = m.L, m.ndims
    return _spatial_profile(Δ, L, ndims)
end

function _spatial_profile(Δ, L, ndims)
    # check to see if we are already in spatial profile form 
    if Base.ndims(Δ) == 3 || Base.ndims(Δ) == 4
        return Δ
    end

    if ndims == 2
        N = L * L
        if length(Δ) == N # swave 
            evs = zeros(5, L, L)
            evs[1, :, :] = reshape(Δ, L, L)
            return evs
        else # dwave 
            evs = zeros(5, L, L)
            for (n, i) in enumerate(1:N:5N)
                evi = Δ[i:(i+N-1)]
                evs[n, :, :] = reshape(evi, L, L)
            end
            return evs
        end

    elseif ndims == 3
        evs = zeros(7, L, L, L)
        N = L * L * L
        if length(Δ) == N # s wave case 
            evs[1, :, :, :] = reshape(Δ, L, L, L)
            return evs
        else
            for (n, i) in enumerate(1:N:7N)
                evi = Δ[i:(i+N-1)]
                evs[n, :, :, :] = reshape(evi, L, L, L)
            end
            return evs
        end
    else
        println("$ndims not implemented yet")
    end
end

function plot_spatial_profile(m::ModelParams; Δ, title=nothing)
    if Base.ndims(Δ) == 3
        return _plot_spatial_profile(m, Δ=Δ, title=title)
    elseif Base.ndims(Δ) == 4
        ps = []
        for slice in 1:L
            pslice = _plot_spatial_profile(m, Δ=Δ, title=title, slice=slice)
            title!(pslice, "z=$slice")
            push!(ps, pslice)
        end
        title = "J=$(m.J), V0=$(m.V0), V1=$(m.V1), θ=$(θ_to_π(m.θ)) on $(m.L)×$(m.L)×$(m.L) lattice"
        p = plot(ps..., layout=Plots.grid(3, 5,
                widths=[1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5]), size=(1500, 1000), plot_title=title)
        return p
    end
end

function _plot_spatial_profile(m::ModelParams; Δ, title=nothing, slice::Int=1)
    evs = spatial_profile(m, Δ=Δ)
    if Base.ndims(Δ) == 4
        evs = evs[:, :, :, slice]
    end

    p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false, aspect_ratio=true)
    for x in 1:L
        for y in 1:L

            # bonds 
            plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[5, x, y]), alpha=10 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)

            # onsite dot 
            scatter!(p, [x], [y], ms=100 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
        end
    end

    xlabel!(p, "Site (x)")
    ylabel!(p, "Site, (y)")
    if !isnothing(title)
        title!(p, title)
    elseif ndims == 2
        title!(p, "J=$(m.J), V0=$(m.V0), V1=$(m.V1), θ=$(θ_to_π(m.θ)) on $(m.L)×$(m.L) lattice")
    end
    return p
end

function plot_in_configuration_space(m::ModelParams; Δ, title::Union{String,Nothing}=nothing, slice::Int=1)

    Δ = spatial_profile(m, Δ=Δ)
    N = numsites(m)

    if m.ndims == 2
        titles = ["on-site", "-x̂ bond", "ŷ bond", "x̂ bond", "-ŷ bond"]
    elseif m.ndims == 3
        titles = ["on-site", "-x̂ bond", "ŷ bond", "x̂ bond", "-ŷ bond",
            "ẑ bond", "-ẑ bond"]
    else
        println("sorry, hun")
        return
    end

    ps = [] # list of plots 
    for b in 1:length(titles)
        p = plot(xlims=(0, 2 * π), ylims=(0, 2 * π), grid=false, aspect_ratio=true)

        if ndims == 2
            Δbond = Δ[b, :, :]

            for r0 in 1:N
                ϕ1, ϕ2 = site_to_configuration_space(r0; m=m)
                x, y = site_to_coordinate(r0, m=m)
                col = Δbond[x, y] < 0 ? "blue" : "red"
                scatter!(p, [ϕ1], [ϕ2], ms=100 * abs(Δbond[x, y]), c=col, legend=:false)
            end

        elseif ndims == 3
            Δbond = Δ[b, :, :, slice]

            for r0 in 1:N
                ϕ1, ϕ2, ϕ3 = site_to_configuration_space(r0; m=m)
                x, y, z = site_to_coordinate(r0, m=m)
                if z == slice
                    col = Δbond[x, y] < 0 ? "blue" : "red"
                    scatter!(p, [ϕ1], [ϕ2], ms=100 * abs(Δbond[x, y]), c=col, legend=:false)
                end
            end

        else
            println("Tough luck")
            return
        end

        xlabel!(p, "ϕ₁")
        ylabel!(p, "ϕ₂")
        title!(p, titles[b])
        push!(ps, p)
    end

    if ndims == 2
        push!(ps, plot_spatial_profile(m, Δ=Δ, title="Real space"))
        p = plot(ps..., layout=Plots.grid(2, 3,
                widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1000), plot_title="Configuration space")
    elseif ndims == 3
        push!(ps, plot_spatial_profile(m, Δ=Δ[:, :, :, slice], title="Real space"))
        p = plot(ps..., layout=Plots.grid(2, 4,
                widths=[1 / 4, 1 / 4, 1 / 4]), size=(1500, 1000), plot_title="Configuration space")
    end

    return p
end

function ΔLGE_to_ΔBdG(m::ModelParams; Δ)

    L, ndims = m.L, m.ndims
    N = numsites(m)

    # get the spatial profile out of the flat thing 
    evs = spatial_profile(m, Δ=Δ)

    Δ_BdG = zeros(N, N)
    for r in 1:N
        coords = site_to_coordinate(r, m=m)

        if ndims == 2
            # coordinates
            x, y = site_to_coordinate(r, m=m)

            # get the nearest neighbours 
            rL, rU, rR, rD = nearest_neighbours(r, m=m)

            # fill in the Δ matrix 
            Δ_BdG[r, r] = evs[1, x, y] # on-site 
            Δ_BdG[r, rL] = evs[2, x, y] # left 
            Δ_BdG[r, rU] = evs[3, x, y] # up 
            Δ_BdG[r, rR] = evs[4, x, y] # right
            Δ_BdG[r, rD] = evs[5, x, y] # down
        elseif ndims == 3
            # coordinates
            x, y, z = site_to_coordinate(r, m=m)

            # get the nearest neighbours 
            rL, rU, rR, rD, rzU, rzD = nearest_neighbours(r, m=m)

            # fill in the Δ matrix 
            Δ_BdG[r, r] = evs[1, x, y, z] # on-site 
            Δ_BdG[r, rL] = evs[2, x, y, z] # left 
            Δ_BdG[r, rU] = evs[3, x, y, z] # up 
            Δ_BdG[r, rR] = evs[4, x, y, z] # right
            Δ_BdG[r, rD] = evs[5, x, y, z] # down
            Δ_BdG[r, rzU] = evs[6, x, y, z] # ẑ up 
            Δ_BdG[r, rzD] = evs[7, x, y, z] # ẑ down
        else
            println("Bummer. $ndims dimensions not yet implemented")
        end
    end

    return Δ_BdG
end

function ΔBdG_to_ΔLGE(m::ModelParams; Δ)
    ndims, L = m.ndims, m.L
    # and now, go backwards 
    if m.ndims == 2
        evs_rec = zeros(5, L, L)
    elseif m.ndims == 3
        evs_rec = zeros(7, L, L, L)
    else
        println("Sorry, $ndims dims not available")
    end
    N = numsites(m)
    for r in 1:N

        if ndims == 2

            # coordinates
            x, y = site_to_coordinate(r, m=m)

            # get the nearest neighbours 
            rL, rU, rR, rD = nearest_neighbours(r, m=m)

            evs_rec[1, x, y] = Δ[r, r] # on-site 
            evs_rec[2, x, y] = Δ[r, rL] # left 
            evs_rec[3, x, y] = Δ[r, rU]# up 
            evs_rec[4, x, y] = Δ[r, rR] # right
            evs_rec[5, x, y] = Δ[r, rD] # down


        elseif ndims == 3

            # coordinates
            x, y, z = site_to_coordinate(r, m=m)

            # get the nearest neighbours 
            rL, rU, rR, rD, rzU, rzD = nearest_neighbours(r, m=m)

            # fill in the Δ matrix 
            evs_rec[1, x, y, z] = Δ[r, r] # on-site 
            evs_rec[2, x, y, z] = Δ[r, rL]# left 
            evs_rec[3, x, y, z] = Δ[r, rU]# up 
            evs_rec[4, x, y, z] = Δ[r, rR] # right
            evs_rec[5, x, y, z] = Δ[r, rD]# down
            evs_rec[6, x, y, z] = Δ[r, rzU] # ẑ up 
            evs_rec[7, x, y, z] = Δ[r, rzD] # ẑ down
        else
            println("No can do, sorry")
            return
        end
    end
    return evs_rec
end

function ΔBdG_to_ΔLGE_flat(m::ModelParams; Δ)
    ndims = m.ndims
    evs_rec = ΔBdG_to_ΔLGE(m, Δ=Δ)
    N = numsites(m)

    if ndims == 2
        evs_flat = zeros(5 * N)
        for (n, i) in enumerate(1:N:(5N))
            evi = evs_rec[n, :, :]
            evs_flat[i:(i+N-1)] = reshape(evi, N)
        end
    elseif ndims == 3
        evs_flat = zeros(7 * N)
        for (n, i) in enumerate(1:N:(7N))
            evi = evs_rec[n, :, :, :]
            evs_flat[i:(i+N-1)] = reshape(evi, N)
        end
    else
        println("Sucks to suck. Sorry :(")
    end

    return evs_flat
end

function to_5N_LGE_Δ(maxev; L)
    @assert length(maxev) == 3 * L * L
    evs = zeros(3, L, L)
    for (n, i) in enumerate(1:(L*L):(3*L*L))
        evi = maxev[i:(i+L*L-1)]
        evs[n, :, :] = reshape(evi, L, L)
    end

    # now put in the missing blocks
    Δx = evs[2, :, :]
    Δy = evs[3, :, :]
    Δxminus = circshift(Δx, (-1, 0)) # mapping between +x̂ and -x̂
    Δyminus = circshift(Δy, (0, 1)) # mapping between +ŷ and -ŷ

    # construct new eigenvector 
    new_evs = zeros(5, L, L)
    new_evs[1, :, :] = evs[1, :, :] # on-site
    new_evs[2, :, :] = Δx
    new_evs[3, :, :] = Δy
    new_evs[4, :, :] = Δxminus
    new_evs[5, :, :] = Δyminus

    evs_flat = zeros(5 * L * L)
    for (n, i) in enumerate(1:(L*L):(5*L*L))
        evi = new_evs[n, :, :]
        evs_flat[i:(i+L*L-1)] = reshape(evi, L * L)
    end

    return evs_flat
end

function to_7N_LGE_Δ(maxev; L)
    N = L * L * L
    @assert length(maxev) == 4 * N
    evs = zeros(4, L, L, L) .* 1im
    for (n, i) in enumerate(1:N:4*N)
        evi = maxev[i:(i+N-1)]
        evs[n, :, :, :] = reshape(evi, L, L, L)
    end

    # now put in the missing blocks
    Δx = evs[2, :, :, :]
    Δy = evs[3, :, :, :]
    Δz = evs[4, :, :, :]
    Δxminus = circshift(Δx, (-1, 0, 0)) # mapping between +x̂ and -x̂
    Δyminus = circshift(Δy, (0, 1, 0)) # mapping between +ŷ and -ŷ
    Δzminus = circshift(Δz, (0, 0, 1)) # mapping between +ẑ and -ẑ

    new_evs = zeros(7, L, L, L) .* 1im
    new_evs[1, :, :, :] = evs[1, :, :, :] # on-site term 
    new_evs[2, :, :, :] = Δx
    new_evs[3, :, :, :] = Δy
    new_evs[4, :, :, :] = Δxminus
    new_evs[5, :, :, :] = Δyminus
    new_evs[6, :, :, :] = Δz
    new_evs[7, :, :, :] = Δzminus

    evs_flat = zeros(7 * N) .* 1im
    for (n, i) in enumerate(1:N:7N)
        evi = new_evs[n, :, :, :]
        evs_flat[i:(i+N-1)] = reshape(evi, N)
    end

    return evs_flat
end

function _symmetry_character(m::ModelParams; Δ)
    L, ndims = m.L, m.ndims

    # rotation about the ẑ axis 
    Δrot = zeros(size(Δ)...)
    for x in 1:L
        for y in 1:L
            for n in 1:5
                x̃ = -y + L + 1
                ỹ = x
                Δrot[n, x̃, ỹ] = Δ[n, x, y]
            end
        end
    end

    # on-site terms
    os = dot(Δ[1, :, :], Δrot[1, :, :])

    # now for the bonds
    bds = dot(Δrot[3, :, :], Δ[2, :, :]) + dot(Δrot[4, :, :], Δ[3, :, :]) + dot(Δrot[5, :, :], Δ[4, :, :]) + dot(Δrot[2, :, :], Δ[5, :, :])

    return os + bds
end

function symmetry_character(m::ModelParams; Δ)
    Δ = spatial_profile(m, Δ=Δ)
    Δ ./= norm(Δ)

    L, ndims = m.L, m.ndims
    @assert ndims == 2 "There is a bug with the 3D sol'n"

    if ndims == 2
        return _symmetry_character(m, Δ=Δ)
    elseif ndims == 3
        @warn "We are assuming a rotation about the ẑ axis"
        χs = []
        for slice in 1:L
            Δslice = Δ[:, :, :, slice]
            push!(χs, _symmetry_character(m, Δ=Δslice))
        end
        return χs
    else
    end
end

# function symmetry_character(m::ModelParams; Δ)
#     L, ndims = m.L, m.ndims
#     if Base.ndims(Δ) == 1
#         Δ = spatial_profile(m, Δ=Δ)
#     end
#     Δ ./= norm(Δ)

#     Δrot = zeros(size(Δ)...)
#     for x in 1:L
#         for y in 1:L
#             for n in 1:5
#                 x̃ = -y + L + 1
#                 ỹ = x
#                 Δrot[n, x̃, ỹ] = Δ[n, x, y]
#             end
#         end
#     end

#     # on-site terms
#     os = dot(Δ[1, :, :], Δrot[1, :, :])

#     # now for the bonds
#     bds = dot(Δrot[3, :, :], Δ[2, :, :]) + dot(Δrot[4, :, :], Δ[3, :, :]) + dot(Δrot[5, :, :], Δ[4, :, :]) + dot(Δrot[2, :, :], Δ[5, :, :])

#     return os + bds
# end

function finite_minimum(matrix)
    fins = matrix[isfinite.(matrix)]
    if prod(size(fins)...) > 0
        return minimum(fins)
    else
        return NaN
    end
end

function finite_maximum(matrix)
    fins = matrix[isfinite.(matrix)]
    if prod(size(fins)...) > 0
        return maximum(fins)
    else
        return NaN
    end
end

function save_structs(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    h5open(path, "w") do file
        for fn in fnames
            n = Name(fn)
            d = getfield(struc, fn)
            write(file, n, d)
        end
    end
end

function load_diagonalized_H(loadpath::String)
    f = h5open(loadpath, "r")
    d = read(f)
    close(f)
    return DiagonalizedHamiltonian(d["L"], d["t"], d["Q"],
        d["μ"], d["θ"], d["ϕx"], d["ϕy"], d["ϕz"], d["J"],
        d["periodic"], d["ndims"], d["disorder"], d["E"], d["U"])
end

function bin_results(arr; nbins)
    δbin = ceil(Int, length(arr) / nbins)
    red = zeros(nbins)
    for (bin, idx) in enumerate(1:δbin:length(arr))
        if idx + δbin > length(arr)
            subarr = arr[idx:end]
        else
            subarr = arr[idx:idx+δbin]
        end
        red[bin] = mean(subarr)
    end
    return red
end

function hist_counts(arr; nbins::Int, normalize::Bool=false)
    δbin = (maximum(arr) - minimum(arr)) / nbins
    counts = zeros(nbins)
    for (bin, idx) in enumerate(minimum(arr)+δbin:δbin:maximum(arr))
        prev = idx - δbin
        subarr = arr[arr.>=prev.&&arr.<idx]
        counts[bin] = length(subarr)
    end
    if normalize
        counts /= δbin
    end
    return counts
end