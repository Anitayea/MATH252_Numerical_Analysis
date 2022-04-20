import BandedMatrices
using LinearAlgebra, Arpack
import SparseArrays
import Plots

#Sets matrix A with the most efficient way, making it into sparsearray to save time and space
function setA(n)
    A = SparseArrays.spzeros(n,n) 
    for i = 1:n
        A[i,i] = 2
        if(i>1)
            A[i,i-1] = -1
        end
        if(i<n)
            A[i,i+1] = -1
        end
    end
    return A
end

#Alternative for setting matrix A, not used 
function setAA(n) #assume n>=3 or else need exception for assigning matrix
    M = SparseArrays.spzeros(n,n)
    for i=1:n
        for j=1:n
            if (j==i)
                M[i,i] = 2.0
            elseif (i!=n && j==i+1)
                M[i,j]= -1.0
            elseif (i!=1 && j==i-1)
                M[i,j]= -1.0
            else
                continue
            end
        end
    end
    return M
end


#setb as a matrix of 1, and multiplies 1/n^2 to move the h to this side to make it easier to implement
setb(n) = (1.0/n^2)*ones(n,1)


function dul(A,n)
    D = SparseArrays.spzeros(size(A))
    U = SparseArrays.spzeros(size(A))
    L = SparseArrays.spzeros(size(A))
    for i in 1:size(A)[1]
        D[i,i] = A[i,i]
        if(i>1)
            L[i,i-1] = A[i,i-1]
        end
        if(i<n)
            U[i,i+1] = A[i,i+1]
        end
    end
    return D,U,L
end

#Alternative for D, not used
function setD(A)
    n = size(A)[1]
    M = SparseArrays.spzeros(n,n)
    for i=1:n
        M[i,i] = A[i,i]
    end
    return M
end

#set M and N from the equation on the textbook
setM(D,L,w) = D/w+L
setN(D,U,w) = (1-w)/w*D - U

function relaxation(n)
    #epsilon for error to stop criterion
    e = 10^-8
    #set up A, b, D, U, L
    A = setA(n)
    b = setb(n)
    #=
    If used the setD, we need this 
    D = setD(A)
    U = triu(A) - D
    L = tril(A) - D
    =#
    D, U, L = dul(A, n)

    #use r to calculate the optimal omega for relaxation method to make it much more efficient
    #
    r = maximum(cos(k*pi/(n+1)) for k in 1:n)
    w = 2/(1+sqrt(1-r^2))

    M = setM(D,L,w)
    N = setN(D,U,w)

    #initialize the x starting with all zeros
    x1 = zeros(n)
    #x2 = zeros(n)
    #use k to keep track of the iteration times
    k=0
    #initialize a vector nR to keep track of the residues
    nR = zeros(0)
    # use the while condition to stop criterion
    while (norm(A*x1-b)/norm(b) >= e)
        #increment k by 1 as iteration times
        k+=1
        #relaxtaion 
        x2 = M\(N * x1 + b) 
        x1 = x2
        #append residues to nR
        append!(nR,norm(A*x1-b))
    end
    return x1,k,nR
end


T = @time relaxation(5000)
#println(T)

# take the variables accordingly and plot the residue and iteration
x1, k, residue = relaxation(5000)
x = 1:1:k
y = residue
Plots.plot(x,y, xlabel = "iterations", ylabel = "residue") 

