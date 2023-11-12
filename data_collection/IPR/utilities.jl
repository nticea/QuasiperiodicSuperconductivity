using DataFrames, CSV
include("../../src/results.jl")

function load_dfs(; dirname="data")
    # read files 
    files = readdir(joinpath(@__DIR__, dirname))
    df = DataFrame(L=[], μ=[], J=[], ϕx=[], ϕy=[], ϕz=[],
        ipr_real=[], ipr_k=[], E=[], pot=[])

    for f in files
        if endswith(f, ".csv")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, dirname, f)))
                dfi = convert_df_arrays(dfi, "ipr_real")
                dfi = convert_df_arrays(dfi, "ipr_k")
                dfi = convert_df_arrays(dfi, "E")
                append!(df, dfi)
            catch e
            end
        end
    end
    dropmissing!(df)

    return df
end