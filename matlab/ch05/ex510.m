% pressdrop.m 
di = 0.75e-2; z = 12; T = 418.15; P1 = 510; mu = 13.8e-6; M = 18; R = 8.314e-3; x0 = 400; 
G = 20:0.1:50; n = length(G); 
for k = 1:n    
Re = di*G(k)/mu; f = 1/(0.79*log(Re) - 0.64)^2;    
fP = @(x) x.^2 + (G(k)^2*R*T/M)*(f*z/di + 2*log(P1./x)) - P1^2; % define eqn.    
P2(k) = fsolve(fP,x0); % use the built-in solver fsolve to solve the eqn. 
end 
plot(G, P2), grid, xlabel('G(kg/sec/m^2)'), ylabel('P_1-P_2 (kPa)')  