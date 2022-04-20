using LinearAlgebra
using Plots

x = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0]
y = [0.6761488864859304; 0.6345697680852508; 0.6396283580587062; 0.6132010027973919;
0.5906142598705267; 0.5718728461471725; 0.5524549902830562; 0.538938885654085;
0.5373495476994958; 0.514904589752926; 0.49243437874655027]


function jacob(x,y)
    ab = [1,1]
    iteration  =  0
    J = zeros(2,2)
    F = zeros(2)
    while iteration<=100
    G_1, G_2, H_11, H_12, H_21, H_22 = 0,0,0,0,0,0
    for i in 1:length(x)
        G_1 += (2*ab[1] - 2*ab[2] * y[i] - 2*x[i]*y[i]) / (ab[2]+x[i])^2
        G_2 += ab[1] * (-2*ab[1] + 2*ab[2]*y[i] + 2*x[i]*y[i]) / (ab[2]+x[i])^3
        H_11 += 2 / (ab[2] +x[i])^2
        H_12 += (-4*ab[1] + 2*ab[2] * y[i] + 2*x[i] * y[i]) / (ab[2] + x[i])^3
        H_21 += (-4*ab[1] + 2*ab[2] * y[i] + 2*x[i] * y[i]) / (ab[2] + x[i])^3
        H_22 += (6*ab[1]^2 - 4*ab[1] * ab[2] * y[i] - 4*ab[1] * x[i] * y[i]) / (ab[2] + x[i])^4
    end
    F[1], F[2] = G_1, G_2
    J[1,1], J[1,2], J[2,1], J[2,2] = H_11, H_12, H_21, H_22
    ab = ab - inv(J)*F
    iteration += 1
    end
    return ab
end

function plot(x)
    M = zeros(0)
    l = []
    for e in x
        f = ab[1]/(ab[2]+e)
        append!(l,f)
    end
    return l
end

ab = jacob(x,y)
l = plot(x)

Plots.plot(x,l)
Plots.scatter!(x,y)