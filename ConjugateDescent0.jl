using LinearAlgebra
import SparseArrays
import Plots
import PlotlyJS

# Domain size
Lx, Ly = 2, 1

# Number of discretization points along x and y, including the boundary points
nx, ny = 101, 101

function discretize(nx, ny)
    hx, hy = Lx/(nx - 1), Ly/(ny - 1)
    Dxx = (1/hx^2) * Tridiagonal(-ones(nx-3), 2ones(nx-2), -ones(nx-3))
    Dyy = (1/hy^2) * Tridiagonal(-ones(ny-3), 2ones(ny-2), -ones(ny-3))
    A = kron(Dxx, I(ny-2)) + kron(I(nx-2), Dyy)
    xgrid = Lx/(nx-1) * (1:nx-2)
    ygrid = Ly/(ny-1) * (1:ny-2)
    x_2d = reshape([x for y in ygrid, x in xgrid], (nx-2)*(ny-2))
    y_2d = reshape([y for y in ygrid, x in xgrid], (nx-2)*(ny-2))
    b = sin.(4π*x_2d) + sin.(2π*y_2d)
    return SparseArrays.SparseMatrixCSC(A), b
end

function plot_solution(f)
    f = reshape(f, ny-2, nx-2)

    # Boundary condition
    z = [zeros(nx)'; zeros(ny-2) f zeros(ny-2); zeros(nx)']
    xgrid = Lx/(nx-1) * (0:nx-1)
    ygrid = Ly/(ny-1) * (0:ny-1)

    # Plots.contourf(xgrid, ygrid, z, c=:viridis, levels=50)
    PlotlyJS.plot(PlotlyJS.contour(x=xgrid, y=ygrid, z=z))
end

# Calculate matrix and right-hand side of linear system
A, b = discretize(nx, ny)

########################
# Your code comes here #
e = 10^(-8)
function CG(A,b)
    x = zeros(size(A)[1])
    # x is a matrix
    r = A*x-b #2n
    # r is a vector
    d = r
    # d is a vector 
    i = 0

    # d,r,x,b vectors
    # A matrix
    # w, beta scalars

    while (norm(r)/norm(b)>=e)
        da = d'*A
        # d'A: vector*matrix: (2a-1)n = 2an-n
        # da is a vector
        dad = da*d
        # da*d: vector*vector: (2n-1)
        # sum: (2a-1)n+(2n-1) = 2an-n+2n-1 = 2an+n-1
        # dad is a scalar
        w = (d'*r)/(dad) 
        # (d'*r): vector*vector: (2n-1)
        # d'*r becomes a scalar
        # (d'*r)/dad: scalar/scalar : 1
        # sum: 2n-1+1 = 2n
        x = x-w*d
        # w*d: scalar*vector: n
        # w*d is a vector
        # x-(w*d): vector - vector : n
        # x is a vector
        # sum: n+n = 2n
        r = A*x-b

        #Ax^{k+1} - b = A(x^{k} - ω_k d_k) - b = r^{k} - ω_k A d_k :)


        # A*x: matrix*vector: (2a-1)n
        # A*x is a vector
        # (A*x)-b: vector - vector: n
        # sum = 2an-n+n = 2an
        beta = (da*r)/(dad)
        # da*r: vector*vector: (2n-1)
        # da*r is a scalar
        # (da*r)/dad: scalar/scalar: 1
        # sum = 2n-1+1 = 2n
        d = r-beta*d
        # beta*d: scalar*vector: n
        # beta*d is a vector
        # r - (beta*d): vector-vector: n
        # sum = 2n
        i+=1
        # i not included
        # sumall = 2an-n+2an+n-1+2n+2n+2an+2n+2n = 6an+8n-1
    end
    return x,i
end

(x,i) = CG(A,b)
f = A\b
display(x)
########################

# Plot the solution
plot_solution(f)
plot_solution(x)