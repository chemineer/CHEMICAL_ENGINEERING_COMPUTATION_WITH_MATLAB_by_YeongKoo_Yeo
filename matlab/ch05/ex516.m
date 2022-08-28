% eqplen.m: calculate the equivalent length of pipe fittings and total length of the pipe 
cfdat = zeros(14,3); cfdat(2,1) = 6; cfdat(11,1) = 2; 
cfdat(:,2) = 100*[10 8 5 3 5 10 15 10 8 15 1.5 1 1.6 0]'; 
cfdat(:,3) = 0.1*[3 2 1.5 1 1.5 2.5 40 20 2.5 15 0.5 0 5 10]'; 
g = 32.2; d = 6.065; D = d/12; pr = 1.5e-4; w = 75000; mu = 1.25; rho = 64.3; % data 
Lst = 78; z = 8; Area = pi*D^2/4; Nre = 6.31*w/(d*mu); v = 0.0509*w/(rho*d^2); Q = w/(8.02*rho); 
if Nre <= 2100 % laminar flow    
f = 64/Nre; % Darcy friction factor 
else % turbulent: friction factor by Chen eqn.    
Av = pr/3.7/D + (6.7/Nre)^0.9; f = 4./(-4*log10(pr/D/3.7 - 5.02*log10(Av)/Nre)).^2; 
end 
Ksum = sum(cfdat(1:12,1).*cfdat(1:12,2)); Klsum = sum(cfdat(1:12,1).*cfdat(1:12,3)); 
Kt = Ksum/Nre + Klsum*(1 + 1/d); 
Kee1 = sum(cfdat(13:14,1).*cfdat(13:14,2)); Kee2 = sum(cfdat(13:14,1).*cfdat(13:14,3)); 
K = Kee1/Nre + Kee2 + Kt; Kt = Kt + K + f*Lst/D; Leq = K*D/f, L = Lst + Leq 
if Nre <= 2100, delP = 0.0034*mu*w/(d^4*rho); 
else, delP = 0.000336*f*w^2/(d^5*rho); end 
delPsi = delP*L/100 + z*rho/144; delH = 0.000483*f*L*w^2/(d^5*rho^2) + z; 