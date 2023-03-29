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

