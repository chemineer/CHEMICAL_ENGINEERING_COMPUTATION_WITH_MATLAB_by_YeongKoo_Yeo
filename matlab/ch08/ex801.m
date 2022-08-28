% ht1D.m: one-dimensional heat transfer (conduction and radiation)
A = 1; T2 = 700; Ta = 1273; sigma = 5.676e-8; xspan = [0 0.2]; T0 = 290; delT = 10;
criT = 1e-3; k = 1;
while delT > criT
[x T] = ode45(@slab1T, xspan, T0, [], A,T2,Ta,sigma);
delT = abs(T2 - T(end)); T2 = T(end); k = k+1;
end
plot(x,T), xlabel('x(m)'), ylabel('T(K)'), grid
k, T(end)
function dT = slab1T(x,T,A,T2,Ta,sigma)
% heat transfer within a one-dimensional slab
dT = - sigma*(T2^4 - Ta^4)/(30*(1 + 0.002*T)*A);
end
