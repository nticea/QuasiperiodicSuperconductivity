## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots, StatsPlots
using CSV, DataFrames, Statistics
using ForwardDiff, Polynomials, LsqFit, Loess


include("../../src/model.jl")
include("../../src/meanfield.jl")
include("../../src/results.jl")
include("utilities.jl")

## PARAMETERS ## 

ndims = 3
Q = (√5 - 1) / 2
θ = π / 7
savefigs = false
figpath = mkpath(joinpath(@__DIR__, "figures"))
T_cutoff = 1e-2

# read files 
files = readdir(joinpath(@__DIR__, "data"))
df = load_dfs()
df = df[(df.θ.==θ).&(df.Q.==Q).&(df.ndims.==ndims), :]

# drop L=27 for now
df = df[(df.L.!=27), :]

# set a χ* so that we can determine Tc 
χstar = 0.25

function Tc_and_deriv(arr, Ts, val)
    Tc = NaN
    try
        idxlist = sortperm(arr)
        x = arr[idxlist]
        y = log10.(Ts[idxlist])
        itpTc = LinearInterpolation(x, y)
        Tclog10 = itpTc(val)
        Tc = 10^Tclog10
    catch e
    end

    # derivative  
    dχdT = NaN
    try
        idxlist = sortperm(Ts)
        y = arr[idxlist]
        x = log10.(Ts[idxlist])
        itpderiv = LinearInterpolation(x, y)
        dχdT = ForwardDiff.derivative(itpderiv, Tclog10)
    catch e
    end

    return Tc, dχdT
end

function LinearSpline(x, y, t)
    Y = Array{Float64}(undef, length(t))
    for i = 1:length(t)
        k = searchsortedfirst(x, t[i])
        if k <= 1
            k = 2
        elseif k >= length(x)
            k = length(x) - 1
        end
        λ = (t[i] - x[k-1]) / (x[k] - x[k-1])
        Y[i] = λ * (y[k] - y[k-1]) + y[k-1]
    end
    return Y
end

function calculate_dχdT(arr, Ts)
    idxlist = sortperm(Ts)
    y = arr[idxlist]
    x = log10.(Ts[idxlist])
    x̂ = LinRange(minimum(x), maximum(x), 1000)

    model = loess(x, y, span=0.5)
    data = Loess.predict(model, x̂)
    window_size = length(x̂ / 20)
    ŷ = [mean(data[i:min(i + window_size - 1, end)]) for i in 1:length(data)]
    return diff(ŷ), x̂[2:end]

    # deg = 5
    # model = Polynomials.fit(x, y, deg)
    # ŷ = model.(x̂)
    # model(x, p) = p[1] ./ (1.0 .+ exp.(-(x .- p[2]) ./ p[3])) .+ p[4]
    # initial_guess = [maximum(y), median(x), 1.0, minimum(y)]
    # fit = LsqFit.curve_fit(model, x, y, initial_guess)
    # params = coef(fit)
    # ŷ = model(x̂, params)

    # @show x
    # @show y
    # @show x̂
    # @show ŷ

    # return diff(ŷ), x̂[2:end]

    # return diff(y), x[2:end]
    # t = LinRange(-0.1, 1.1, 6000)
    # F(t) = LinearSpline(x, y, t)
    # ŷ = F(x̂)

    # return diff(ŷ), x̂[2:end]

    # @show x
    # @show y
    # @assert 1 == 0
    # itpderiv = LinearInterpolation(x, y)
    # xnew = LinRange(minimum(x), maximum(x), 300)
    # dχdT = []
    # xs = []
    # for xn in xnew
    #     try
    #         dχ = ForwardDiff.derivative(itpderiv, xn)
    #         push!(dχdT, dχ)
    #         push!(xs, xn)
    #     catch
    #     end
    # end
    # return dχdT, xs
end

Ls = unique(df.L)
df_Tc = DataFrame(L=[], J=[], θ=[], Q=[], ndims=[],
    Tc_swave=[], Tc_dwave=[], dχ_swave=[], dχ_dwave=[],
    dχ_swave_all=[], dχ_dwave_all=[], Ts_swave=[], Ts_dwave=[])
for L in Ls
    dfL = df[(df.L.==L), :]
    Js = unique(dfL.J)
    for J in Js
        dfJ = dfL[(dfL.J.==J), :]
        χs, Ts = dfJ.χ, dfJ.T
        χswave, χdwave, Ts = get_χswave_dwave(χs, Ts)

        Tc_swave, dχ_swave = Tc_and_deriv(χswave, Ts, χstar)
        Tc_dwave, dχ_dwave = Tc_and_deriv(χdwave, Ts, χstar)

        dχ_swave_all, swave_Ts = calculate_dχdT(χswave, Ts)
        dχ_dwave_all, dwave_Ts = calculate_dχdT(χdwave, Ts)

        dfi = DataFrame(L=[L], J=[J], θ=[θ], Q=[Q], ndims=[ndims],
            Tc_swave=[Tc_swave], Tc_dwave=[Tc_dwave],
            dχ_swave=[dχ_swave], dχ_dwave=[dχ_dwave],
            dχ_swave_all=[dχ_swave_all], dχ_dwave_all=[dχ_dwave_all],
            Ts_swave=[swave_Ts], Ts_dwave=[dwave_Ts])
        append!(df_Tc, dfi)
    end
end

# plotting 
Tc_p, dχ_p = plot(xlabel="J", ylabel="Tc"), plot(xlabel="J", ylabel="dχ/dlog10T")
cmapb = cgrad(:blues, length(Ls), categorical=true)
cmapr = cgrad(:reds, length(Ls), categorical=true)

for (l, L) in enumerate(Ls)
    dfL = df_Tc[(df_Tc.L.==L), :]
    Js, Tcs, Tcd, dχs, dχd = dfL.J, dfL.Tc_swave, dfL.Tc_dwave, dfL.dχ_swave, dfL.dχ_dwave

    plot!(Tc_p, Js, Tcs, c=cmapr[l], label=nothing)
    scatter!(Tc_p, Js, Tcs, c=cmapr[l], label="s-wave, L=$L")
    plot!(Tc_p, Js, Tcd, c=cmapb[l], label=nothing)
    scatter!(Tc_p, Js, Tcd, c=cmapb[l], label="d-wave, L=$L")

    plot!(dχ_p, Js, dχs, c=cmapr[l], label=nothing)
    scatter!(dχ_p, Js, dχs, c=cmapr[l], label="s-wave, L=$L")
    plot!(dχ_p, Js, dχd, c=cmapb[l], label=nothing)
    scatter!(dχ_p, Js, dχd, c=cmapb[l], label="d-wave, L=$L")
end

ps = []
pswave = plot(title="s-wave", xlabel="log10T", ylabel="dχ/dlog10T")
pdwave = plot(title="d-wave", xlabel="log10T", ylabel="dχ/dlog10T")
Js = unique(df.J)
cmap = cgrad(:viridis, length(Js), categorical=true)

for (j, J) in enumerate(Js)
    p = plot(xlabel="T", ylabel="dχ/dlog10T")
    for (l, L) in enumerate(Ls)
        dfL = df_Tc[(df_Tc.L.==L).&(df_Tc.J.==J), :]
        if size(dfL)[1] > 0
            dχ_swave_all, dχ_dwave_all, Ts_swave, Ts_dwave = dfL.dχ_swave_all, dfL.dχ_dwave_all, dfL.Ts_swave, dfL.Ts_dwave
            plot!(p, Ts_swave, dχ_swave_all, c=cmapr[l], label=nothing)
            #scatter!(p, Ts_swave, dχ_swave_all, c=cmapr[l], label="s-wave, L=$L")
            plot!(p, Ts_dwave, dχ_dwave_all, c=cmapb[l], label=nothing)
            #scatter!(p, Ts_dwave, dχ_dwave_all, c=cmapb[l], label="d-wave, L=$L")

            plot!(pswave, Ts_swave, dχ_swave_all, c=cmap[j], label="J=$J")
            plot!(pdwave, Ts_dwave, dχ_dwave_all, c=cmap[j], label=nothing)
            #hmap = heatmap(zeros(1, 1), xlims=xlims(pdwave), ylims=ylims(pdwave), clims=(minimum(Js), maximum(Js)), cmap=:viridis, alpha=0)
        end
    end
    push!(ps, p)
end

# make a colourbar 
p2 = plot(pswave, pdwave, layout=Plots.grid(1, 2,
        widths=[1 / 2, 1 / 2]), size=(1700, 800), plot_title="dχ/dlog10T with Q=$(round(Q, digits=3)), θ=$(θ_to_π(θ))")

# plot to make: T*(dchi/dT) as a function of T, plotted for all the different L on same plot 
# make a different plot for each J 
