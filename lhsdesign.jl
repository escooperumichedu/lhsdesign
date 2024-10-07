function lhsdesign(n::Int, p::Int)
    # returns a Latin hypercube sample matrix of size n-by-p. 
    # For each column of X, the n values are randomly distributed with one from each interval (0,1/n), (1/n,2/n), ..., (1 - 1/n,1), and randomly permuted.

    # Empty parameter space
    X = zeros(n, p)
    
    # n - samples
    # p - variables

    # Window intervals
    a = (0:1:n-1) / n # Left hand side
    b = (1:1:n) / n # Right hand side


    for j in 1:p
        for i in 1:n
            X[i, j] = a[i] + (b[i] - a[i]) * rand()
        end
        X[:, j] .= shuffle!(X[:, j])
    end

    return X
end
