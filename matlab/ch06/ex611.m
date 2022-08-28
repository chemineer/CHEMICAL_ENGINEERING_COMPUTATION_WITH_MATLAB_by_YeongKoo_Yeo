% exocstr.m 
% Data: 
F0 = 40; F = 40; Fj = 49.9;Ca0=0.55; V = 48; rho = 50; rhoj = 62.3; Cp = 0.75; Cj = 1; A = 250; U = 150;  
T0 = 530; Tj0 = 530; alp=7.08e10; lam = -3e4; E = 3e4; R = 1.9872; 
% Define nonlinear equations 
fun = @(x) [F0*Ca0 - F*x(1) - alp*V*x(1)*exp(-E/R/x(2));        
rho*Cp*(F0*T0-F*x(2))-lam*alp*V*x(1)*exp(-E/R/x(2))-U*A*(x(2) - x(3));        
rhoj*Cj*Fj*(Tj0 - x(3)) + U*A*(x(2) - x(3))]; 
% Solution of nonlinear equation system 
x0 = [Ca0 T0 Tj0]; % initial guess 
x = fsolve(fun, x0); Ca = x(1), T = x(2), Tj = x(3) 