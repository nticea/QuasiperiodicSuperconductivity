using DataFrames, CSV
include("../../src/results.jl")

function already_computed(df; T, L, Q, θ, ϕx, ϕy, ϕz, J, ndims, Λ)
    # extract only the parameters we are interested in 
    df = df[(df.L.==L).&(df.T.==T).&(df.θ.==θ).&(df.Q.==Q).&(df.ϕx.==ϕx).&(df.ϕy.==ϕy).&(df.ϕz.==ϕz).&(df.ndims.==ndims).&(df.J.==J).&(df.Λ.==Λ), :]
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
        T=[], Λ=[], χ=[], dχdlogT=[])
    for f in files
        if endswith(f, ".csv")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
                dfi = convert_df_arrays(dfi, "χ")
                dfi = convert_df_arrays(dfi, "dχdlogT")
                append!(df, dfi)
            catch e
            end
        end
    end
    dropmissing!(df)

    return df
end

function get_χswave_dwave(χs, Ts)
    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χs = χs[sortidx]
    χs = [reshape(χ, 4, 4) for χ in χs]

    if length(Ts) > 0
        # on-site
        χswave = [χ[1, 1] for χ in χs]

        # make the d-wave components 
        xx, yy = [χ[2, 2] for χ in χs], [χ[3, 3] for χ in χs]
        xy, yx = [χ[2, 3] for χ in χs], [χ[3, 2] for χ in χs]
        if ndims == 2
            χdwave = xx + yy - xy - yx
        elseif ndims == 3
            zz = [χ[4, 4] for χ in χs]
            xz, zx = [χ[2, 4] for χ in χs], [χ[4, 2] for χ in χs]
            yz, zy = [χ[3, 4] for χ in χs], [χ[4, 3] for χ in χs]
            χdwave = xx + yy + zz - xy - yx - xz - zx - yz - zy
        else
            println("sorry")
            χdwave = nothing
        end
    end

    return χswave, χdwave, Ts
end

function get_χswave(χs, Ts)
    sortidx = sortperm(Ts)
    Ts = Ts[sortidx]
    χs = χs[sortidx]
    χs = [reshape(χ, 4, 4) for χ in χs]

    if length(Ts) > 0
        χswave = [χ[1, 1] for χ in χs]
    end

    return χswave, Ts
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