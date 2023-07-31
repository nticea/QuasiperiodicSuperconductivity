using DataFrames, CSV
include("../../src/results.jl")

function already_computed(df; L, Q, θ, ϕx, ϕy, J)
    # extract only the parameters we are interested in 
    df = df[(df.L.==L).&(df.θ.==θ).&(df.Q.==Q).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy).&(df.J.==J), :]
    if size(df)[1] > 0
        println("Already computed!")
        return true
    end
    return false
end

function load_dfs()
    # read files 
    files = readdir(joinpath(@__DIR__, "data"))
    df = DataFrame(L=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[],
        T=[], χ=[])
    for f in files
        if endswith(f, ".csv")
            dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
            dfi = convert_df_arrays(dfi, "χ")
            append!(df, dfi)
        end
    end

    return df
end