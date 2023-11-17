using DataFrames, CSV
include("../../src/results.jl")

function already_computed(df, m::ModelParams; λ, E₀)
    # extract only the parameters we are interested in 
    df = df[(df.L.==m.L).&(df.θ.==m.θ).&(df.Q.==m.Q).&(df.ϕx.==m.ϕx).&(df.ϕy.==m.ϕy).&(df.ϕz.==m.ϕz).&(df.ndims.==m.ndims).&(df.disorder.==m.disorder).&(df.J.==m.J).&(df.λ.==λ).&(df.E₀.==E₀), :]
    if size(df)[1] > 0
        println("Already computed!")
        return true
    end
    return false
end

function save_results(m::ModelParams; λ, E₀, α₀, ipr_real, ipr_k)
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH:MM:SS")
    datapath = joinpath(@__DIR__, "data")
    mkpath(datapath)
    savepath = joinpath(datapath, "$(L)L_$(J)J" * timestamp * ".csv")
    df = DataFrame(L=[m.L], J=[m.J], Q=[m.Q], θ=[m.θ],
        ϕx=[m.ϕx], ϕy=[m.ϕy], ϕz=[m.ϕz], ndims=[m.ndims],
        λ=[λ], E₀=[E₀], α₀=[α₀], ipr_real=[ipr_real],
        ipr_k=[ipr_k], disorder=[m.disorder])
    CSV.write(savepath, df)
    flush(stdout)
end

function load_dfs()
    # read files 
    files = readdir(joinpath(@__DIR__, "data"))
    df = DataFrame(L=[], J=[], Q=[], θ=[], ϕx=[], ϕy=[], ϕz=[],
        ndims=[], disorder=[], λ=[], E₀=[], α₀=[], ipr_real=[], ipr_k=[])
    for f in files
        if endswith(f, ".csv")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
                append!(df, dfi)
            catch
            end
        end
    end
    dropmissing!(df)
    return df
end