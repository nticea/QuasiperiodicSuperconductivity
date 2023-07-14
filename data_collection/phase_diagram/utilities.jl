function already_computed(; L, t, Q, μ, θ, ϕx, ϕy, V0, V1, J, periodic)
    # read files 
    files = readdir(joinpath(@__DIR__, "data"))
    df_LGE = DataFrame(L=Int64[], J=Float64[], Q=Float64[], θ=Float64[],
        ϕx=Float64[], ϕy=Float64[], V0=Float64[], V1=Float64[],
        T=Float64[], λ=Float64[], Δ=[], K=[], Π=[])
    for f in files
        if endswith(f, ".csv") && contains(f, "LGE")
            dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
            # process all of the arrays 
            dfi = convert_df_arrays(dfi, "Δ")
            dfi = convert_df_arrays(dfi, "K")
            dfi = convert_df_arrays(dfi, "Π")
            append!(df_LGE, dfi)
        end
    end

    # extract only the parameters we are interested in 
    df_LGE = df_LGE[(df_LGE.L.==L).&(df_LGE.θ.==θ).&(df_LGE.Q.==Q).&(df_LGE.ϕx.==ϕx).&(df_LGE.ϕy.==ϕy).&(df_LGE.V0.==V0).&(df_LGE.V1.==V1).&(df_LGE.J.==J), :]
    if size(df_LGE)[1] > 0
        println("Already computed!")
        return true
    end

    return false
end