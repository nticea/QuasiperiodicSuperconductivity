using DataFrames
include("IndexMatrix.jl")

function update_results!(df::DataFrame; L, λ, J, θ, ϕx, ϕy, V0, V1, T, Δ)
    if isnothing(θ)
        θ = 0
    end
    df2 = DataFrame(L=[L], λ=[λ], J=[J], θ=[θ], ϕx=[ϕx], ϕy=[ϕy], V0=[V0], V1=[V1], T=[T], Δ=[Δ])
    append!(df, df2)
end

function update_results!(df::DataFrame; L, λ, J, Q, θ, ϕx, ϕy, V0, V1, T, Δ, K=zeros(4), Π=zeros(4))
    if isnothing(θ)
        θ = 0
    end
    df2 = DataFrame(L=[L], λ=[λ], J=[J], Q=[Q], θ=[θ], ϕx=[ϕx], ϕy=[ϕy], V0=[V0], V1=[V1], T=[T], Δ=[Δ], K=[K], Π=[Π])
    append!(df, df2)
end

function already_calculated(df::DataFrame; L, J, θ, ϕx, ϕy, V0, V1, T)
    sub = df[(df.L.==L).&(df.J.==J).&(df.V0.==V0).&(df.V1.==V1).&(df.θ.==θ).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy).&(df.T.==T), :]
    return size(sub)[1] > 0
end

function already_calculated(df::DataFrame; L, J, Q, θ, ϕx, ϕy, V0, V1, T)
    sub = df[(df.L.==L).&(df.J.==J).&(df.Q.==Q).&(df.V0.==V0).&(df.V1.==V1).&(df.θ.==θ).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy).&(df.T.==T), :]
    return size(sub)[1] > 0
end

# same as above, but without T (this is for finding Tc using LGE method)
function already_calculated(df::DataFrame; L, J, θ, ϕx, ϕy, V0, V1)
    sub = df[(df.L.==L).&(df.J.==J).&(df.V0.==V0).&(df.V1.==V1).&(df.θ.==θ).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy), :]
    return size(sub)[1] > 0
end

# same as above, but without T (this is for finding Tc using LGE method)
function already_calculated(df::DataFrame; L, J, Q, θ, ϕx, ϕy, V0, V1)
    sub = df[(df.L.==L).&(df.J.==J).&(df.Q.==Q).&(df.V0.==V0).&(df.V1.==V1).&(df.θ.==θ).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy), :]
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

function load_dataframe(path)
    try # try loading the DataFrames
        df = DataFrame(CSV.File(path))

        all_evs_new = []
        for (i, evs) in enumerate(df.Δ)
            # Split the string by multiple delimiters
            result = split(chop(evs; head=1, tail=1), r"[,; ]")
            new_array = result[result.!=""]
            evs_new = parse.(Float64, new_array)
            push!(all_evs_new, evs_new)
        end
        dfcut = df[:, collect(1:size(df)[2]-1)]
        dfcut.Δ = all_evs_new

        return dfcut
    catch error_reading_dataframe # if the file does not exist, create a new dataframe
        @show error_reading_dataframe
        return DataFrame(L=Int64[], J=Float64[], Q=Float64[], θ=Float64[],
            ϕx=Float64[], ϕy=Float64[], V0=Float64[], V1=Float64[],
            T=Float64[], λ=Float64[], Δ=[], K=[], Π=[])
    end
end

function find_Tc_LGE(df::DataFrame; interp_value::Real=1)

    # Get the unique Js, V1s from the dataframe 
    Js, V0s, V1s = unique(df.J), unique(df.V0), unique(df.V1)
    L = df.L[1]

    # Make the new dataframe 
    nodenames = ["L", "J", "V0", "V1", "Tc"]
    Tc_df = DataFrame([name => [] for name in nodenames])

    # For each unique (J,V0) pair... extract the corresponding data across all Ts 
    for J in Js
        for V0 in V0s
            for V1 in V1s

                # Extract the corresponding data across all Ts
                dfsub = df[(df.J.==J).&(df.V0.==V0).&(df.V1.==V1), :]
                λs, Ts = dfsub.λ, dfsub.T

                # Compute the interpolated Tc for this (J,V0) pair
                idxlist = sortperm(λs)
                knots = λs[idxlist]
                if length(knots) > 1
                    Interpolations.deduplicate_knots!(knots, move_knots=true)
                    try
                        interp_linear = linear_interpolation(knots, Ts[idxlist])
                        Tc = interp_linear(interp_value)
                        # Put it into a new dataframe indexed by (J,V0,Tc)
                        df2 = DataFrame(L=[L], Tc=[Tc], J=[J], V0=[V0], V1=[V1])
                        append!(Tc_df, df2)
                    catch e
                        Tc = NaN
                        # Put it into a new dataframe indexed by (J,V0,Tc)
                        df2 = DataFrame(L=[L], Tc=[Tc], J=[J], V0=[V0], V1=[V1])
                        append!(Tc_df, df2)
                    end
                end
            end
        end
    end
    return Tc_df
end

function find_Tc_BdG(df::DataFrame; interp_value::Real=1)

    # Get the unique Js, V1s from the dataframe 
    Js, V0s, V1s = unique(df.J), unique(df.V0), unique(df.V1)
    L = df.L[1]

    # Make the new dataframe 
    nodenames = ["L", "J", "V0", "V1", "Tc"]
    Tc_df = DataFrame([name => [] for name in nodenames])

    # For each unique (J,V0) pair... extract the corresponding data across all Ts 
    for J in Js
        for V0 in V0s
            for V1 in V1s

                # Extract the corresponding data across all Ts
                dfsub = df[(df.J.==J).&(df.V0.==V0).&(df.V1.==V1), :]
                Δs, Ts = maximum.(dfsub.Δ), dfsub.T

                # Compute the interpolated Tc for this (J,V0) pair
                idxlist = sortperm(Δs)
                knots = Δs[idxlist]
                if length(knots) > 1
                    Interpolations.deduplicate_knots!(knots, move_knots=true)
                    try
                        interp_linear = linear_interpolation(knots, Ts[idxlist])
                        Tc = interp_linear(interp_value)
                        # Put it into a new dataframe indexed by (J,V0,Tc)
                        df2 = DataFrame(L=[L], Tc=[Tc], J=[J], V0=[V0], V1=[V1])
                        append!(Tc_df, df2)
                    catch e
                        Tc = NaN
                        # Put it into a new dataframe indexed by (J,V0,Tc)
                        df2 = DataFrame(L=[L], Tc=[Tc], J=[J], V0=[V0], V1=[V1])
                        append!(Tc_df, df2)
                    end
                end
            end
        end
    end
    return Tc_df
end

function θ_to_π(θ)
    if isnothing(θ)
        return "Untilted"
    end
    for n in 1:10
        if isapprox(n * θ, π)
            return "π/$n"
        end
    end
    return θ
end

function colour_phase(x1::Int, x2::Int, x3::Int; all_evs, numpts::Int=10)
    if ndims(all_evs) == 3
        val = all_evs[x1, x2, x3]
    elseif ndims(all_evs) == 2
        val = all_evs[x2, x3]
    end
    if val < 0
        return "blue"
    else
        return "red"
    end
end

function spatial_profile(maxev; L)
    evs = zeros(5, L, L)
    for (n, i) in enumerate(1:(L*L):(5*L*L))
        evi = maxev[i:(i+L*L-1)]
        evs[n, :, :] = reshape(evi, L, L)
    end
    return evs
end


function plot_spatial_profile(evs; L, title=nothing)

    if ndims(evs) != 3
        evs = spatial_profile(evs, L=L)
    end

    p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
    for x in 1:L
        for y in 1:L

            # bonds 
            plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
            plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)

            # onsite dot 
            scatter!(p, [x], [y], ms=100 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)
        end
    end

    xlabel!(p, "Site (x)")
    ylabel!(p, "Site, (y)")
    if !isnothing(title)
        title!(p, title)
    end
    return p
end

function plot_in_config_space(toplot; L::Int, Q::Real, θ::Union{Real,Nothing}, title::Union{String,Nothing})

    @assert ndims(toplot) == 2

    p = plot(xlims=(0, 2 * π), ylims=(0, 2 * π), grid=false, aspect_ratio=true)
    for x in 1:L
        for y in 1:L
            ϕ1, ϕ2 = coordinate_to_configuration_space(x, y, L=L, Q=Q, θ=θ)
            if toplot[x, y] < 0
                col = "blue"
            else
                col = "red"
            end

            scatter!(p, [ϕ1], [ϕ2], ms=100 * abs(toplot[x, y]), c=col, legend=:false)
        end
    end

    xlabel!(p, "ϕ₁")
    ylabel!(p, "ϕ₂")
    if !isnothing(title)
        title!(p, title)
    end
    return p
end

function plot_LGE_Δ(df; idx)
    L = df.L[idx]
    J = df.J[idx]
    V0 = df.V0[idx]
    V1 = df.V1[idx]
    T = df.T[idx]
    θ = θ_to_π(df.θ[idx])
    ϕx = θ_to_π(df.ϕx[idx])
    ϕy = θ_to_π(df.ϕy[idx])

    if length(df.Δ[1]) == 5 * L^2
        symmetry = "d-wave"
    elseif length(df.Δ[1]) == L^2
        symmetry = "s-wave"
    else
        @error "Symmetry not recognized"
        @assert 1 == 0
    end

    maxev = df.Δ[idx]
    λ = df.λ[idx]
    if symmetry == "d-wave"
        evs = zeros(5, L, L)
        for (n, i) in enumerate(1:(L*L):(5*L*L))
            evi = maxev[i:(i+L*L-1)]
            evs[n, :, :] = reshape(evi, L, L)
        end
    elseif symmetry == "s-wave"
        evs = reshape(maxev, L, L)
    end

    p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
    if λ > 0 && symmetry == "d-wave"
        for x in 1:L
            for y in 1:L

                # bonds 
                plot!(p, [x, x - 1], [y, y], lw=10 * abs(evs[1, x, y]), alpha=10 * abs(evs[1, x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
                plot!(p, [x, x], [y, y + 1], lw=10 * abs(evs[2, x, y]), alpha=10 * abs(evs[2, x, y]), c=colour_phase(2, x, y, all_evs=evs), legend=:false)
                plot!(p, [x, x + 1], [y, y], lw=10 * abs(evs[3, x, y]), alpha=10 * abs(evs[3, x, y]), c=colour_phase(3, x, y, all_evs=evs), legend=:false)
                plot!(p, [x, x], [y, y - 1], lw=10 * abs(evs[4, x, y]), alpha=10 * abs(evs[4, x, y]), c=colour_phase(4, x, y, all_evs=evs), legend=:false)

                # onsite dot 
                if abs.(maximum(evs[5, x, y])) > 1e-6
                    scatter!(p, [x], [y], ms=100 * abs(evs[5, x, y]), c=colour_phase(5, x, y, all_evs=evs), legend=:false)
                end

            end
        end
    elseif λ > 0 && symmetry == "s-wave"
        for x in 1:L
            for y in 1:L
                if abs.(maximum(evs[x, y])) > 1e-6
                    scatter!(p, [x], [y], ms=100 * abs(evs[x, y]), c=colour_phase(1, x, y, all_evs=evs), legend=:false)
                end
            end
        end
    end

    xlabel!(p, "Site (x)")
    ylabel!(p, "Site, (y)")
    title!(p, "T=$(round(T,digits=4)), λ=$(round(λ,digits=2)) \n Δ(J=$J, θ=$θ, ϕx=$ϕx, ϕy=$ϕy, \n V0=$V0, V1=$(round(V1,digits=2)))", fontsize=4)
    return p
end

function plot_all_Δs(loadpath)
    df = load_dataframe(loadpath)
    df = sort(df, :V0)

    J = df.J[1]
    V1 = df.V1[1]

    # for every unique V1, find the Δ with the λ closest to 0 
    V0s = unique(df.V0)
    global hmaps = []
    for V0 in V0s
        # get the corresponding data
        subdf = df[(df.V0.==V0), :]
        λs = subdf.λ
        idx = argmin(abs.(λs .- 1))
        Ts = subdf.T

        @show J, V0, V1
        @show Ts
        @show λs

        println("")
        push!(hmaps, plot_LGE_Δ(subdf; idx=idx))
    end
    p = plot(hmaps..., layout=Plots.grid(3, 3, widths=[1 / 3, 1 / 3, 1 / 3]), size=(1500, 1500), aspect_ratio=:equal)
end

function ΔLGE_to_ΔBdG(Δ_LGE; L::Int)
    evs = zeros(5, L, L)
    for (n, i) in enumerate(1:(L*L):(5*L*L))
        evi = Δ_LGE[i:(i+L*L-1)]
        evs[n, :, :] = reshape(evi, L, L)
    end

    Δ_BdG = zeros(L * L, L * L)
    for x in 1:L
        for y in 1:L
            # coordinate 
            r = coordinate_to_site(x, y, L=L)

            # get the nearest neighbours 
            rL, rU, rR, rD = nearest_neighbours(r, L=L)

            # fill in the Δ matrix 
            Δ_BdG[r, rL] = evs[1, x, y] # left 
            Δ_BdG[r, rU] = evs[2, x, y] # up 
            Δ_BdG[r, rR] = evs[3, x, y] # right
            Δ_BdG[r, rD] = evs[4, x, y] # down
            Δ_BdG[r, r] = evs[5, x, y] # on-site 
        end
    end

    return Δ_BdG
end

function ΔBdG_to_ΔLGE(Δ_BdG; L::Int)
    # and now, go backwards 
    evs_rec = zeros(5, L, L)

    for r in 1:(L*L)
        # coordinates
        x, y = site_to_coordinate(r, L=L)

        # get the nearest neighbours 
        rL, rU, rR, rD = nearest_neighbours(r, L=L)

        evs_rec[1, x, y] = Δ_BdG[r, rL] # left 
        evs_rec[2, x, y] = Δ_BdG[r, rU]# up 
        evs_rec[3, x, y] = Δ_BdG[r, rR] # right
        evs_rec[4, x, y] = Δ_BdG[r, rD] # down
        evs_rec[5, x, y] = Δ_BdG[r, r] # on-site 
    end

    return evs_rec
end

function ΔBdG_to_ΔLGE_flat(Δ_BdG; L::Int)
    evs_rec = ΔBdG_to_ΔLGE(Δ_BdG, L=L)

    evs_flat = zeros(5 * L * L)
    for (n, i) in enumerate(1:(L*L):(5*L*L))
        evi = evs_rec[n, :, :]
        evs_flat[i:(i+L*L-1)] = reshape(evi, L * L)
    end

    return evs_flat
end

function symmetry_character(Δ; L::Int)
    if ndims(Δ) == 1
        Δ = spatial_profile(Δ, L=L)
    end

    # normalize 
    Δ ./= norm(Δ)

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
    os = dot(Δ[5, :, :], Δrot[5, :, :])

    # now for the bonds
    bds = dot(Δrot[2, :, :], Δ[1, :, :]) + dot(Δrot[3, :, :], Δ[2, :, :]) + dot(Δrot[4, :, :], Δ[3, :, :]) + dot(Δrot[1, :, :], Δ[4, :, :])

    return os + bds

end

finite_minimum(matrix) = minimum(matrix[isfinite.(matrix)])
finite_maximum(matrix) = maximum(matrix[isfinite.(matrix)])
