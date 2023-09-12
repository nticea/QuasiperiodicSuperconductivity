
# Define a function to split a matrix into LxL squares
function split_into_squares(matrix, L)
    # Get the dimensions of the input matrix
    rows, cols = size(matrix)

    # Check if the matrix can be evenly divided into squares
    if rows % L != 0 || cols % L != 0
        throw(ArgumentError("Matrix dimensions are not divisible by L"))
    end

    # Initialize an array to store the squares
    squares = Matrix{Matrix}(undef, (Int(rows ÷ L), Int(cols ÷ L)))

    # Loop through the matrix and extract each square
    for i in 1:Int(rows ÷ L)
        for j in 1:Int(cols ÷ L)
            row_start = Int((i - 1) * L + 1)
            row_end = Int(i * L)
            col_start = Int((j - 1) * L + 1)
            col_end = Int(j * L)
            squares[i, j] = matrix[row_start:row_end, col_start:col_end]
        end
    end

    return squares
end

# Define a function to split a matrix into LxLxL cubes
function split_into_cubes(matrix, L)
    # Get the dimensions of the input matrix
    rows, cols, height = size(matrix)

    # Check if the matrix can be evenly divided into squares
    if rows % L != 0 || cols % L != 0 || height % L != 0
        throw(ArgumentError("Matrix dimensions are not divisible by L"))
    end

    # Initialize an array to store the squares
    cubes = Array{Any}(undef, (Int(rows ÷ L), Int(cols ÷ L), Int(height ÷ L)))

    # Loop through the matrix and extract each square
    for i in 1:Int(rows ÷ L)
        for j in 1:Int(cols ÷ L)
            for k in 1:Int(height ÷ L)
                row_start = Int((i - 1) * L + 1)
                row_end = Int(i * L)
                col_start = Int((j - 1) * L + 1)
                col_end = Int(j * L)
                height_start = Int((k - 1) * L + 1)
                height_end = Int(k * L)
                cubes[i, j, k] = matrix[row_start:row_end, col_start:col_end, height_start:height_end]
            end
        end
    end

    return cubes
end
