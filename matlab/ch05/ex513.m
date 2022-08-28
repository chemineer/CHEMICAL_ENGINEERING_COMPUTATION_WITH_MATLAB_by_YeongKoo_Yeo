% sol3tank.m 
clear all; 
% Data 
L = 100*[6 2 1]; d = [0.1524 0.0625 0.1016]; 
ep = 5e-6*ones(1,3); Z = [8 35 50]; rho = 1e3; mu = 1e-3; g = 9.8; Q1 = 0.042; 
% x(1)=Ws; x(2)=v2, x(3)=v3, x(4)=f1, x(5)=f2, x(6)=f3 
x0 = [10 5 5 1e-3 1e-3 1e-3]; % initial condition 
x = fsolve(@pipsys,x0,[],L,d,ep,Z,rho,mu,g,Q1); % Solve nonlinear equation system 
Ws = x(1); v2 = x(2); v3 = x(3); Q2 = pi*v2*d(2)^2 / 4; Q3 = pi*v3*d(3)^2 / 4; 
fprintf('Power of the pump = %g m of column of water\n', Ws); 
fprintf('Volumetric flow rate through pipe 1 = %g m^3/sec\n', Q1) 
fprintf('Volumetric flow rate through pipe 2 = %g m^3/sec\n', Q2) 
fprintf('Volumetric flow rate through pipe 3 = %g m^3/sec\n', Q3) 