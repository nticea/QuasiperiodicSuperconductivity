import Base.setindex!
import Base.getindex

struct IndexMatrix
    data::Matrix{}
    rows
    cols
end

function Base.setindex!(s::IndexMatrix, value, row::Real, col::Real)
    cidx = argmin(abs.(col .- s.cols))
    ridx = argmin(abs.(row .- s.rows))
    s.data[ridx, cidx] = value
end

function Base.getindex(s::IndexMatrix, row::Real, col::Real)
    cidx = argmin(abs.(col .- s.cols))
    ridx = argmin(abs.(row .- s.rows))
    return s.data[ridx, cidx]
end

function IndexMatrix(rowvals, colvals)
    function make_axes(arr)
        sort!(arr)
        maxval, minval = maximum(arr), minimum(arr)
        s = [arr[n] - arr[n-1] for n in 2:length(arr)]
        spacing = minimum(s)
        newvals = collect(minval:spacing:maxval)
        return newvals
    end
    cols, rows = make_axes(colvals), make_axes(rowvals)
    data = fill(NaN, (length(rows), length(cols)))
    return IndexMatrix(data, rows, cols)
end