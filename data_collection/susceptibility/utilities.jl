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

function fit_χ(Ts, χs; nbits::Real=3)
    # split up Ts into 5 bits 
    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χs = χs[sortidx]

    wl = ceil(Int, length(Ts) / nbits) # window length 
    errs = []
    for idx in 1:(length(Ts)-wl)
        T = log10.(Ts[idx:idx+wl])
        χ = χs[idx:idx+wl]
        b, a = linear_fit(T, χ)
        χ̂ = a .* T .+ b
        err = norm(χ - χ̂)
        push!(errs, err)
    end

    minerr = argmin(errs)
    T = log10.(Ts[minerr:minerr+wl])
    χ = χs[minerr:minerr+wl]
    b, a = linear_fit(T, χ)
    χ̂ = a .* T .+ b

    return 10 .^ T, χ̂, a, b, errs[minerr]
end