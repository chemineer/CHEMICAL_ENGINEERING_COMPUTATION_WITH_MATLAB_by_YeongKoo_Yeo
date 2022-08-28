% sufrad.m
radt.A1 = 10; radt.A2 = 15; radt.W = 1; radt.H = 1; radt.L = 10;
radt.sig = 5.67e-8; radt.T1 = 1000; radt.T3 = 300; radt.ep1 = 0.9;
radt.ep2 = 0.5; radt.q2 = 77100;
% x(1) = J1, x(2) = J2, x(3) = T2
x0 = [500 500 500];  % initial guess
x = fsolve(@frad, x0, [], radt)
q1 = radt.ep1*radt.A1*(radt.sig*(radt.T1)^4 - x(1))/(1-radt.ep1)
q2 = -radt.ep2*radt.A2*(radt.sig*x(3)^4 - x(2))/(1-radt.ep2)

function f = frad(x,radt) 
% x(1) = J1, x(2) = J2, x(3) = T2 
A1 = radt.A1; A2 = radt.A2; W = radt.W; H = radt.H; L = radt.L; sig = radt.sig; 
T1 = radt.T1; T3 = radt.T3; ep1 = radt.ep1; ep2 = radt.ep2; q2 = radt.q2; 
X = W/H; Y = L/H;  
F12 = (2/(pi*X*Y))*(log(sqrt((1+X^2)*(1+X^2)/(1+X^2+Y^2))) +... 
X*sqrt(1+Y^2)*atan(X/sqrt(1+Y^2)) + Y*sqrt(1+X^2)*atan(Y/sqrt(1+X^2)) - X*atan(X) - Y*atan(Y)); 
F13 = 1 - F12; F23 = F13*A1/A2;  
f = [ep1*A1*(sig*T1^4 - x(1))/(1-ep1) - F12*A1*(x(1)-x(2)) - F13*A1*(x(1)- sig*T3^4); 
ep2*A2*(sig*x(3)^4 - x(2))/(1-ep2) - F12*A1*(x(2)-x(1)) - F23*A2*(x(2)-sig*T3^4);    
q2 + ep2*A2*(sig*x(3)^4 - x(2))/(1-ep2)]; 
end 