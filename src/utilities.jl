
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

function uniform_susceptibility_components(χ; ndims::Int)
    if Base.ndims(χ) == 1
        if ndims == 2
            χ = reshape(χ, 3, 3)
        elseif ndims == 3
            χ = reshape(χ, 4, 4)
        else
            println("do you know who i am")
        end
    end
    # on-site
    χswave = χ[1, 1]

    # make the d-wave components 
    xx, yy = χ[2, 2], χ[3, 3]
    xy, yx = χ[2, 3], χ[3, 2]
    if ndims == 2
        χdwave = xx + yy - xy - yx
    elseif ndims == 3
        zz = χ[4, 4]
        xz, zx = χ[2, 4], χ[4, 2]
        yz, zy = χ[3, 4], χ[4, 3]
        χdwave = xx + yy + zz - xy - yx - xz - zx - yz - zy
    else
        println("sorry")
        χdwave = nothing
    end

    return χswave, χdwave
end