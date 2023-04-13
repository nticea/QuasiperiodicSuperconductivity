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

function update_results!(df::DataFrame; L, λ, J, V0, T)
    df2 = DataFrame(L=[L], λ=[λ], J=[J], V0=[V0], T=[T])
    append!(df, df2)
end

function already_calculated(df::DataFrame; L, J, V0, T)
    sub = df[(df.L.==L).&(df.J.==J).&(df.V0.==V0).&(df.T.==T), :]
    return size(sub)[1] > 0
end

function load_dataframe(path)
    try # try loading the DataFrames
        return DataFrame(CSV.File(path))
    catch error_reading_dataframe # if the file does not exist, create a new dataframe
        @show error_reading_dataframe
        nodenames = ["L", "J", "V0", "T", "λ"]
        return DataFrame([name => [] for name in nodenames])
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

    # Get the unique Js, V0s from the dataframe 
    Js, V0s = unique(df.J), unique(df.V0)
    L = df.L[1]

    # Make the new dataframe 
    nodenames = ["L", "J", "V0", "Tc"]
    Tc_df = DataFrame([name => [] for name in nodenames])

    # For each unique (J,V0) pair... extract the corresponding data across all Ts 
    for J in Js
        for V0 in V0s

            # Extract the corresponding data across all Ts
            dfsub = df[(df.J.==J).&(df.V0.==V0), :]
            λs, Ts = dfsub.λ, dfsub.T

            # Compute the interpolated Tc for this (J,V0) pair
            knots = reverse(λs)
            if length(knots) > 1
                Interpolations.deduplicate_knots!(knots, move_knots=true)
                try
                    interp_linear = linear_interpolation(knots, reverse(Ts))
                    Tc = interp_linear(interp_value)
                    # Put it into a new dataframe indexed by (J,V0,Tc)
                    df2 = DataFrame(L=[L], Tc=[Tc], J=[J], V0=[V0])
                    append!(Tc_df, df2)
                catch e
                    Tc = NaN
                    # Put it into a new dataframe indexed by (J,V0,Tc)
                    df2 = DataFrame(L=[L], Tc=[Tc], J=[J], V0=[V0])
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

