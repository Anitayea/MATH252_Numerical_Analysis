#1.6
#using recurion NOT array
function b(n)
    if(n==0)
        return 0
    else
        return n%2+10*b(n÷2) 
    end
end
#using arrat
function to_binary(n)
    A = zeros(Int, 0)
    while n>0
        append!(A, n%2)
        n ÷= 2
    end
    return A
end
println(to_binary(10))

#bitstring: used to check which to take x+ or x-