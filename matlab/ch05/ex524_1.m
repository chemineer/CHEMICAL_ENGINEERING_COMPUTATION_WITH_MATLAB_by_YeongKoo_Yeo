% dpgasflowc.m: pressure drop for critical flow 
clear all; 
G = 582/3600; d1 = 0.10226; d2 = 0.01; p1 = 1.2e6; Fa = 1; rho = 10.25; mu = 1.3e-5; 
gam = 1.3; K = 1; 
beta = d2/d1; A1 = pi*d1^2/4; q1 = G/rho; v1 = q1/A1; % upstream velocity (m/s) 
Ao = pi*d2^2/4; % orifice area 
Nre = d1*v1*rho/mu; % Reynolds number 
Cd = 0.73; % critical flow 
Y = G/(Cd*Fa*Ao)/sqrt(2*rho*p1*gam*(2/(1+gam))^((gam+1)/(gam-1))); 
dp = p1/(0.9953 + 0.9054/sqrt(K) + 0.1173/K - 0.0195/K^1.5); 
fprintf('Expansion factor = %g\n', Y); fprintf('Pressure drop = %g kPa\n', dp/1000) 