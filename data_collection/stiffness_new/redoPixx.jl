## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using CSV
using DataFrames
using StatsPlots

## PARAMETERS ## 

L = 17 # the full system is L × L 
Q = (√5 - 1) / 2
θ = π / 7
V0 = 1
V1 = -1.5

files = readdir(joinpath(@__DIR__, "data"))
df_BdG = DataFrame(L=Int64[], J=Float64[], Q=Float64[], θ=Float64[],
    ϕx=Float64[], ϕy=Float64[], V0=Float64[], V1=Float64[],
    T=Float64[], λ=Float64[], Δ=[], K=[], Π=[])
df_LGE = copy(df_BdG)
for f in files
    if endswith(f, ".csv")
        dfi = DataFrame(CSV.File(joinpath(@__DIR__, "data", f)))
        # process all of the arrays 
        dfi = convert_df_arrays(dfi, "Δ")
        dfi = convert_df_arrays(dfi, "K")
        dfi = convert_df_arrays(dfi, "Π")
        if contains(f, "BdG")
            append!(df_BdG, dfi)
        elseif contains(f, "LGE")
            append!(df_LGE, dfi)
        end
    end
end

# Open the file in read mode
file = open("example.txt", "r")

# Read the entire contents of the file
contents = read(file, String)

# Close the file
close(file)

# Define the string to search for
search_string = "pattern"

# Use regular expressions to find matches
matches = matchall(r"\b$search_string\b", contents)

# Print the matches
for match in matches
    println(match.match)
end