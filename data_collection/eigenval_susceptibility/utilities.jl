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
        T=[], χswave=[], χdwave=[], dχswave=[], dχdwave=[], Δswave=[], Δdwave=[])#,
    #IPR_real=[], IPR_momentum=[])
    for f in files
        if endswith(f, ".csv")
            try
                dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
                dfi = convert_df_arrays(dfi, "Δswave")
                dfi = convert_df_arrays(dfi, "Δdwave")
                #dfi = get_IPRs(df)
                append!(df, dfi)
            catch e
                println("file $f could not be loaded...")
                println(e)
            end
        end
    end
    dropmissing!(df)

    return df
end

function get_IPRs!(df::DataFrame)

    all_IPRS_real = []
    all_IPRS_momentum = []
    for row in eachrow(df)
        L, ndims = row.L, row.ndims

        IPRS_real = [IPR_real(row.Δswave)]
        IPRS_momentum = [IPR_momentum(row.Δswave, ndims, L)]

        Δ = reshape_eigenvector(row.Δdwave, ndims, L)
        Δ = _spatial_profile(Δ, L, ndims)

        n = size(Δ)[1]
        for i in 2:n
            if Base.ndims(Δ) == 4
                Δcomp = reshape(Δ[i, :, :, :], prod(size(Δ[i, :, :, :])))
            elseif Base.ndims(Δ) == 3
                Δcomp = reshape(Δ[i, :, :], prod(size(Δ[i, :, :])))
            else
                println("so sad")
            end
            Δcomp ./= norm(Δcomp)
            push!(IPRS_real, IPR_real(Δcomp))
            push!(IPRS_momentum, IPR_momentum(Δcomp, ndims, L))
        end
        push!(all_IPRS_real, IPRS_real)
        push!(all_IPRS_momentum, IPRS_momentum)
    end
    # df.IPR_real = all_IPRS_real
    # df.IPR_momentum = all_IPRS_momentum

    df.IPR_swave_real = [x[1] for x in all_IPRS_real]
    df.IPR_swave_momentum = [x[1] for x in all_IPRS_momentum]
    df.IPR_x_real = [x[2] for x in all_IPRS_real]
    df.IPR_x_momentum = [x[2] for x in all_IPRS_momentum]
    df.IPR_y_real = [x[3] for x in all_IPRS_real]
    df.IPR_y_momentum = [x[3] for x in all_IPRS_momentum]
    df.IPR_z_real = [x[6] for x in all_IPRS_real]
    df.IPR_z_momentum = [x[6] for x in all_IPRS_momentum]

    return df
end

function reshape_eigenvector(vec, ndims, L)
    if ndims == 3
        return to_7N_χ(vec, L)
    elseif ndims == 2
        return to_5N_χ(vec, L)
    else
        println("Worst. day. ever.")
    end
end

function to_5N_χ(maxev, L)
    @assert length(maxev) == 2 * L * L
    evs = zeros(2, L, L)
    for (n, i) in enumerate(1:(L*L):(2*L*L))
        evi = maxev[i:(i+L*L-1)]
        evs[n, :, :] = reshape(evi, L, L)
    end

    # now put in the missing blocks
    Δx = evs[1, :, :]
    Δy = evs[2, :, :]
    Δxminus = circshift(Δx, (-1, 0)) # mapping between +x̂ and -x̂
    Δyminus = circshift(Δy, (0, 1)) # mapping between +ŷ and -ŷ

    # construct new eigenvector 
    new_evs = zeros(5, L, L)
    new_evs[1, :, :] = zeros(size(Δx)) # on-site
    new_evs[2, :, :] = Δx
    new_evs[3, :, :] = Δy
    new_evs[4, :, :] = Δxminus
    new_evs[5, :, :] = Δyminus

    evs_flat = zeros(5 * L * L)
    for (n, i) in enumerate(1:(L*L):(5*L*L))
        evi = new_evs[n, :, :]
        evs_flat[i:(i+L*L-1)] = reshape(evi, L * L)
    end

    return evs_flat
end

function to_7N_χ(maxev, L)
    N = L * L * L
    @assert length(maxev) == 3 * N
    evs = zeros(3, L, L, L)
    for (n, i) in enumerate(1:N:3*N)
        evi = maxev[i:(i+N-1)]
        evs[n, :, :, :] = reshape(evi, L, L, L)
    end

    # now put in the missing blocks
    Δx = evs[1, :, :, :]
    Δy = evs[2, :, :, :]
    Δz = evs[3, :, :, :]
    Δxminus = circshift(Δx, (-1, 0, 0)) # mapping between +x̂ and -x̂
    Δyminus = circshift(Δy, (0, 1, 0)) # mapping between +ŷ and -ŷ
    Δzminus = circshift(Δz, (0, 0, 1)) # mapping between +ẑ and -ẑ

    new_evs = zeros(7, L, L, L)
    new_evs[1, :, :, :] = zeros(size(Δx)) # on-site term 
    new_evs[2, :, :, :] = Δx
    new_evs[3, :, :, :] = Δy
    new_evs[4, :, :, :] = Δxminus
    new_evs[5, :, :, :] = Δyminus
    new_evs[6, :, :, :] = Δz
    new_evs[7, :, :, :] = Δzminus

    evs_flat = zeros(7 * N)
    for (n, i) in enumerate(1:N:7N)
        evi = new_evs[n, :, :, :]
        evs_flat[i:(i+N-1)] = reshape(evi, N)
    end

    return evs_flat
end


# We need to transform each of the eigenvectors into 2D space! 
function fourier_transform(u, ndims::Int, L::Int; minus=false)
    # first, map back to 2D space 
    if ndims == 2
        N = L * L
        u = reshape(u, L, L)
    elseif ndims == 3
        N = L * L * L
        u = reshape(u, L, L, L)
    else
        println("$dims dimensions not implemented")
        return
    end

    # take FT along all spatial dimensions 
    if !minus # this is the expression for U_{q}
        uq = fft(u)
    else # this is the expression for U_{-q}
        uq = conj.(fft(conj.(u)))
    end

    # reshape it back and normalize. FFTW does not normalize!!
    uq = reshape(uq, N) ./ √N

    return uq
end

function IPR_momentum(u, ndims, L)
    # transform to momentum space
    uq = fourier_transform(u, ndims, L)
    u2 = uq .* conj.(uq)
    u4 = u2 .^ 2
    return real(sum(u4))
end

function IPR_real(u)
    u2 = u .* conj.(u)
    u4 = u2 .^ 2
    return real(sum(u4))
end