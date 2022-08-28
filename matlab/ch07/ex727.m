% filtconsP.m: constant pressure filtration
clear all;
t = [4.4, 9.5, 16.3, 24.6, 34.7, 46.1, 59.0, 73.6, 89.4, 107.3]'; % t(s)
V = 1e-3*[0.498, 1.0, 1.501, 2.0, 2.498, 3.002, 3.506, 4.004, 4.502, 5.009]'; % V(m^3)
dP = 338e3; A = 0.0439; cs = 23.47; mu = 8.937e-4;
n = length(t); Y = t./V; C = [V/2 ones(n,1)]; X = inv(C'*C)*C'*Y;
Kp = X(1); B = X(2); alpa = Kp*A^2*dP/(mu*cs); Rm = A*dP*B/mu;
fprintf('Alpha = %g m/kg\n', alpa); fprintf('Rm = %g m^(-1)\n', Rm);
