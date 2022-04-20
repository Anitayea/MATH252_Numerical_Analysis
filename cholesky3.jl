import BandedMatrices
import LinearAlgebra

#=
Computation cost of the unbanded cholesky decomposition is 1/3 n^3, 
as in the code there are 3 loops, ijk, and since j =1:i and k is 1:j-1, it is 1/3*;
while for Lu decomposition it is 2/3 n^3, which is slightly faster. 
To optimize the cholesky decomposition with banded matrices, 
we skip the loops when it exits the band since elements outside of the bands are all Zeros,
and therefore the run time can be lessened to 1/3 k^2 n, 
which is  more efficient for banded matrices, k < n 
=#

function cholesky(M)
    m, n = size(M)
    m!=n && error("Matrix must be square")
    B = BandedMatrices.BandedMatrix(M)
    B.u != B.l && error("Matrix must be symmetric")
    
    n = size(M)[1]
    N = BandedMatrices.BandedMatrix(BandedMatrices.Zeros(n,n),(B.l, B.u))
    b = M.u
    for i=1:n
        for j=1:i
            sum = 0.0
            #skip loop early to optimize the time when the matrix is banded
            if (j<i-b)
                continue
            end
            #
            for k=1:j-1
                #skip loop early to optimize the time when the matrix is banded
                if (k<i-b)
                    continue
                end
                #
                sum += (N[i,k]*N[j,k])
            end
            if(i==j)
                N[i,i] = sqrt((M[i,i])-sum)
            else
                N[i,j] = 1.0 / N[j,j]*(M[i,j]-sum)
            end
        end
        if (N[i,i]<=0)
            println("Not positive definite!")
            return -1
        end
    end
    return N
end

n, u, l = 20000,2,2
A = BandedMatrices.brand(n, u, l)
A = A*A'
T = @time cholesky(A)
println(LinearAlgebra.norm(T*T' - A, Inf))
display(cholesky(A))

B = [6 15 55; 15 55 225; 55 225 979]
B = BandedMatrices.BandedMatrix(B)
t = @time cholesky(B)
println(LinearAlgebra.norm(t*t' - B, Inf))
display(cholesky(B))