% dpgasflow.m: pressure drop for gas flow 
G = 3900/3600; d1 = 0.10226; d2 = 0.035; p1 = 1.2e6; Fa = 1; rho = 10.25; mu = 1.3e-5; gam = 1.3; % data 
beta = d2/d1; A1 = pi*d1^2/4; q1 = G/rho; v1 = q1/A1; % upstream velocity (m/s) 
Nre = d1*v1*rho/mu; % Reynolds number 
Cd = 0.6274 - 0.2354*beta + 0.7858*beta^2; % restriction orifice 
Ao = pi*d2^2/4; % orifice area 
fun = @(x) 1-(0.41+0.35*beta^4)*x/gam/p1 - sqrt(rho*(1-beta^4)/2/x)*G/ (rho*Ao*Cd*Fa); 
x0 = 100; dp = fsolve(fun,x0); 
Y = 1-(0.41+0.35*beta^4)*dp/gam/p1; 
fprintf('Expansion factor = %g\n', Y); fprintf('Pressure drop = %g kPa\n', dp/1000)  