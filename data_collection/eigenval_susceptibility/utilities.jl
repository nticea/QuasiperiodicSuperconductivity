using DataFrames, CSV
include("../../src/results.jl")

function already_computed(df; T, L, Q, θ, ϕx, ϕy, ϕz, J, ndims)
    # extract only the parameters we are interested in 
    df = df[(df.L.==L).&(df.T.==T).&(df.θ.==θ).&(df.Q.==Q).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy).&(df.ϕz.==ϕz).&(df.ndims.==ndims).&(df.J.==J), :]
    if size(df)[1] > 0
        println("Already computed!")
        return true
    end
    return false
end

function load_dfs()
    # read files 
    files = readdir(joinpath(@__DIR__, "data"))
    df = DataFrame(L=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[], ϕz=[], ndims=[],
        T=[], Λ=[], χswave=[], χdwave=[])
    for f in files
        if endswith(f, ".csv")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
                append!(df, dfi)
            catch e
            end
        end
    end
    dropmissing!(df)

    return df
end