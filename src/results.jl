using HDF5

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