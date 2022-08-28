function fv = vwfun(v,T,L,D,rf,dz,dP) 
g = 32.174; gc = 32.174; D = D/12; dP = 144*dP; % psi->lbf/ft^2 
eD = rf./D; % roughness factor 
% density(rho; lbm/ft^3) and viscosity(mu; lbm/ft/s) at T 
rho = 62.122+0.0122*T-(1.54e-4)*T.^2+(2.65e-7)*T.^3-(2.24e-10)*T.^4; 
mu = exp(-11.0318 + 1057.51./(T+214.624)); 
Nre = D.*v.*rho./mu; % Reynolds number 
if Nre < 2100, f = 16./Nre; 
else, den = 16*(log10(eD/3.7 - 5.02*log10(eD/3.7+14.5./Nre)./Nre)).^2; f = 1./ den; end 
fv = v - sqrt((g*dz + gc*dP./rho)./(0.5 - 2*f.*L./D));  
end 