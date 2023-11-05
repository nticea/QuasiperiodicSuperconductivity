using DataFrames, CSV
include("../../src/results.jl")

function already_computed(m::ModelParams, df::DataFrame; T)
    # extract only the parameters we are interested in 
    df = df[(df.L.==m.L).&(df.T.==T).&(df.θ.==m.θ).&(df.Q.==m.Q).&(df.t.==m.t).&(df.μ.==m.μ).&(df.V0.==m.V0).&(df.V1.==m.V1).&(df.ϕx.==m.ϕx).&(df.ϕy.==m.ϕy).&(df.ϕz.==m.ϕz).&(df.ndims.==m.ndims).&(df.J.==m.J), :]
    if size(df)[1] > 0
        for K in df.K
            if any(abs.(K) .> 0)
                println("Already computed!")
                return true
            end
        end
    end
    return false
end

function load_dfs(dirname)
    # read files 
    files = readdir(joinpath(@__DIR__, dirname))
    df = DataFrame(L=[], λ=[], Δ=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[],
        ϕz=[], ndims=[], V0=[], V1=[], t=[], μ=[], periodic=[], T=[],
        K=[], Π=[])
    for f in files
        if endswith(f, ".csv")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, dirname, f)))
                dfi = convert_df_arrays(dfi, "K")
                dfi = convert_df_arrays(dfi, "Π")
                dfi = convert_df_arrays(dfi, "Δ")
                append!(df, dfi)
            catch e
            end
        end
    end
    dropmissing!(df)

    return df
end

function load_LGE(dirname)
    # read files 
    files = readdir(joinpath(@__DIR__, dirname))
    df = DataFrame(L=[], λ=[], Δ=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[],
        ϕz=[], ndims=[], V0=[], V1=[], t=[], μ=[], periodic=[], T=[],
        K=[], Π=[])
    for f in files
        if endswith(f, ".csv") && contains(f, "LGE")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, dirname, f)))
                dfi = convert_df_arrays(dfi, "K")
                dfi = convert_df_arrays(dfi, "Π")
                dfi = convert_df_arrays(dfi, "Δ")
                append!(df, dfi)
            catch e
                @show e
                #println("file $f not valid...?")
            end
        end
    end
    dropmissing!(df)

    return df
end

function load_BdG(dirname)
    # read files 
    files = readdir(joinpath(@__DIR__, dirname))
    df = DataFrame(L=[], λ=[], Δ=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[],
        ϕz=[], ndims=[], V0=[], V1=[], t=[], μ=[], periodic=[], T=[],
        K=[], Π=[])
    for f in files
        if endswith(f, ".csv") && contains(f, "BdG")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, dirname, f)))
                dfi = convert_df_arrays(dfi, "K")
                dfi = convert_df_arrays(dfi, "Π")
                dfi = convert_df_arrays(dfi, "Δ")
                append!(df, dfi)
            catch e
                println("file $f not valid...?")
            end
        end
    end
    dropmissing!(df)

    return df
end