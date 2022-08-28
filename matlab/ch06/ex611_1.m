% multst.m 
% Data: 
F0 = 40; F = 40; Fj = 49.9;Ca0=0.55; V = 48; rho = 50; rhoj = 62.3; Cp = 0.75; Cj = 1; A = 250; U = 150;  
T0 = 530; Tj0 = 530; alp=7.08e10; lam = -3e4; E = 3e4; R = 1.9872; 
% Define nonlinear functions 
fT = @(T) [rho*Cp*(F0*T0 - F*T) - (F0*Ca0*V*lam*alp*exp(-E/R./T))./...    
(F+alp*V*exp(-E/R./T)) - U*A*rhoj*Cj*Fj*(T-Tj0)./(rhoj*Cj*Fj + U*A)]; 
% Solve nonlinear equations by using the built-in solver fzero 
T0 = T0+150; % initial guess 
x = fzero(fT, T0) 
% Plot of f(T) versus T 
Tv = 500:0.1:700; Fv = fT(Tv); Fv0 = Fv*0; plot(Tv,Fv,Tv,Fv0), xlabel('T(R)'), ylabel('f(T)') 