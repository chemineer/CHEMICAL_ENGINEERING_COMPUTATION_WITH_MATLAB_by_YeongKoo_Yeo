function fD = dwfun(D,W,rho,mu,rf,dP,L,K) 
eD = rf./D; % roughness factor 
v = 4*W./(pi*rho*D.^2); % m/s 
Nre = D.*v.*rho./mu; % Reynolds number 
if Nre < 2100    
f = 16./Nre; 
else % Shacham eqn    
den = 16*(log10(eD/3.7 - 5.02*log10(eD/3.7+14.5./Nre)./Nre)).^2; f = 1./den; 
end 
Le = K*D/(4*f); fD = D - ((2*f*(L + Le))./(rho*dP) * (4*W/pi)^2).^0.2;  
end 