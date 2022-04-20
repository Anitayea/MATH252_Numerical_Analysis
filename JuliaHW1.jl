#task 4
function euler_constant(n)
    a=0
    for i = 1:n
        a += 1/i
    end
    return -log(ℯ,n)+ a
end
println(euler_constant(1000))
#or
euler_constant(n) = (-log(ℯ,n)+sum(i=1/i for i in 1:n))
println(euler_constant(1000))

#task 5
function T(n,r)
    if (n<=0||r<=0)
        return 0
    elseif(n==1&&r>=1)
        return 1
    elseif(r==3)
        return 2^n-1
    elseif(r>3&&n>0)
        a=typemax(Int64)
        b=0
        for k = 1:n-1
            b = (2*T(k, r))+T(n-k, r-1)
            if(b<a)
                a=b
            end
        end
        return a
    else
        #when illegal input
        return -1
    end
end

#another method

function T(n,r)
    if (n<=0||r<=0)
        return 0
    elseif (n==1&&r>=1)
        return 1
    elseif (r==3)
        return 2^n-1
    elseif (r>3&&n>0)
        M = [2*T(k, r)+T(n-k, r-1) for k in 1:n-1]
        return minimum(M)
    else
        #when illegal input
        return -1
    end
end

println(T(4,4))
println(T(5,4))
println(T(2,4))
println(T(4,5))
println(T(4,2))