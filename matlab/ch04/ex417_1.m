% pvplot.m 
A = [13.8183 13.8587]; B = [2477.07 2991.32]; C = [233.21 216.64]; 
T = 60; z = 0.6; v = 0:0.01:1; n = length(v); x0 = [0.1 0.6 50]; 
for k = 1:n, x = fsolve(@binflash,x0,[],v(k),T,z,A,B,C);  P(k) = x(3); end 
plot(v,P), xlabel('v(vaporized frac.)'), ylabel('P(kPa)'), grid 