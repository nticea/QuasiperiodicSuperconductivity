## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra

function diis(iterate_func, init_guess, maxiter, conv_tol, diis_size)
    # Define variables
    num_iter = 0
    error_vec = zeros(maxiter)
    f_vec = zeros(diis_size, diis_size)
    b_vec = zeros(diis_size)
    x_vec = zeros(diis_size)

    # Initialize variables
    x = init_guess
    error = Inf
    converged = false

    # Iterate
    while (num_iter < maxiter) && !converged
        # Calculate next iterate
        x_new = iterate_func(x)

        # Calculate error
        error_new = norm(x_new - x)

        # Add new error to vector
        error_vec[num_iter+1] = error_new

        # Check for convergence
        if error_new < conv_tol
            converged = true
        end

        # Add new f vector to matrix
        if num_iter > diis_size
            f_new = zeros(diis_size, diis_size)
            for i = 1:diis_size
                for j = 1:diis_size
                    f_new[i, j] = dot(error_vec[num_iter-diis_size+i:num_iter-1], error_vec[num_iter-diis_size+j:num_iter-1])
                end
            end
            f_vec[1:diis_size-1, 1:diis_size-1] = f_vec[2:diis_size, 2:diis_size]
            f_vec[1:diis_size-1, diis_size] = f_vec[diis_size, 1:diis_size-1]
            f_vec[diis_size, :] = f_new
        end

        # Solve linear system to obtain new x vector
        if num_iter >= diis_size
            b_vec[1:diis_size-1] = b_vec[2:diis_size]
            b_vec[diis_size] = -1
            x_vec = f_vec \ b_vec
            x_new = zero(x)
            for i = 1:diis_size
                x_new += x_vec[i] * iterate_func(error_vec[num_iter-diis_size+i:num_iter-1])
            end
        end

        # Update iterate and error
        x = x_new
        error = error_new

        # Increment iteration count
        num_iter += 1
    end

    # Return final iterate and error
    return x, error
end

iterate_func(x) = 0.5 * (x + 1 / x)

init_guess = 2.0
maxiter = 50
conv_tol = 1e-8
diis_size = 5

x, error = diis(iterate_func, init_guess, maxiter, conv_tol, diis_size)
