using HDF5
using DataFrames

struct Results
    L
    λs
    Js
    V0s
    Ts
end

function save_structs(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    for fn in fnames
        n = Name(fn)
        d = getfield(struc, fn)

        # If the file already exists, then we either append to it or overwrite 
        if isfile(path)
            h5open(path, "r+") do file
                if haskey(file, n) #if key already exists, we want to rewrite 
                    delete_object(file, n)
                    write(file, n, d)
                else
                    write(file, n, d)
                end
            end
        else # If the file does not exist, create it 
            h5open(path, "w") do file
                write(file, n, d)
            end
        end
    end
end

function load_results(loadpath::String)
    f = h5open(loadpath, "r")
    d = read(f)
    return Results(d["L"], d["λs"], d["Js"], d["V0s"], d["Ts"])
end

# function update_results!(df::DataFrame; L, λ, J, V0, T)
#     df2 = DataFrame(L=[L], λ=[λ], J=[J], V0=[V0], T=[T])
#     append!(df, df2)
# end

function update_results!(df::DataFrame; L, λ, J, θ, V0, V1, T, Δ)
    if isnothing(θ)
        θ = 0
    end
    df2 = DataFrame(L=[L], λ=[λ], J=[J], θ=[θ], V0=[V0], V1=[V1], T=[T], Δ=[Δ])
    append!(df, df2)
end

function already_calculated(df::DataFrame; L, J, θ, V0, V1, T)
    sub = df[(df.L.==L).&(df.J.==J).&(df.V0.==V0).&(df.V1.==V1).&(df.θ.==θ).&(df.T.==T), :]
    return size(sub)[1] > 0
end

function load_dataframe(path)
    try # try loading the DataFrames
        df = DataFrame(CSV.File(path))

        all_evs_new = []
        for (i, evs) in enumerate(df.Δ)
            # Split the string by multiple delimiters
            result = split(chop(evs; head=1, tail=1), r",")
            new_array = result[result.!=""]
            evs_new = parse.(Float64, new_array)
            push!(all_evs_new, evs_new)
        end
        dfcut = df[:, collect(1:size(df)[2]-1)]
        dfcut.Δ = all_evs_new

        return dfcut
    catch error_reading_dataframe # if the file does not exist, create a new dataframe
        @show error_reading_dataframe
        nodenames = ["L", "J", "θ", "V0", "V1", "T", "λ", "Δ"]
        # return DataFrame([name => [] for name in nodenames])
        return DataFrame(L=Int64[], J=Float64[], θ=Float64[], V0=Float64[], V1=Float64[],
            T=Float64[], λ=Float64[], Δ=[])
    end
end

function find_Tc(results::Results; interp_value::Real=1)
    L, λs, Js, V0s, Ts = results.L, results.λs, results.Js, results.V0s, results.Ts

    # do interpolations 
    Tcs = zeros(length(Js), length(V0s))
    for k in 1:length(Js)
        for j in 1:length(V0s)
            knots = reverse(λs[k, j, :])
            Interpolations.deduplicate_knots!(knots, move_knots=true)
            try
                interp_linear = linear_interpolation(knots, reverse(Ts))
                Tcs[k, j] = interp_linear(interp_value)
            catch e
                Tcs[k, j] = NaN
            end
        end
    end

    return Tcs
end

function find_Tc(df::DataFrame; interp_value::Real=1)

    # Get the unique Js, V1s from the dataframe 
    Js, V1s = unique(df.J), unique(df.V1)
    L = df.L[1]
    V0 = df.V0[1]

    # Make the new dataframe 
    nodenames = ["L", "J", "V0", "V1", "Tc"]
    Tc_df = DataFrame([name => [] for name in nodenames])

    # For each unique (J,V0) pair... extract the corresponding data across all Ts 
    for J in Js
        for V1 in V1s

            # Extract the corresponding data across all Ts
            dfsub = df[(df.J.==J).&(df.V1.==V1), :]
            λs, Ts = dfsub.λ, dfsub.T

            # Compute the interpolated Tc for this (J,V0) pair
            knots = reverse(λs)
            if length(knots) > 1
                Interpolations.deduplicate_knots!(knots, move_knots=true)
                try
                    interp_linear = linear_interpolation(knots, reverse(Ts))
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
    return Tc_df
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

function plot_LGE_Δ(df; idx)
    L = df.L[idx]
    J = df.J[idx]
    V0 = df.V0[idx]
    V1 = df.V1[idx]
    T = df.T[idx]
    θ = θ_to_π(df.θ[idx])

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

    function colour_phase(x1::Int, x2::Int, x3::Int; all_evs, numpts::Int=10)
        cm = palette([:blue, :red], 2 * numpts + 1)
        val = all_evs[x1, x2, x3]
        max = maximum(abs.(all_evs))
        idx = floor(Int, val / max * numpts + numpts + 1)
        return cm[idx]
    end

    p = plot(xlims=(0, L + 1), ylims=(0, L + 1), grid=false)
    # p = plot(xlims=(0, L + 1), ylims=(-L - 1, 0))
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
    xlabel!(p, "Site (x)")
    ylabel!(p, "Site, (y)")
    title!(p, "T=$T, λ=$(round(abs(λ),digits=2)) \n Δ(J=$J, θ=$θ, V0=$V0, V1=$(round(V1,digits=2)))", fontsize=4)
    return p
end

