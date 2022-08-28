function [phig f] = phigas(state,eos,T,P,Tc,Pc,w) 
% Estimates the fugacity coefficient using the virial and cubic EOS. 
% input 
% state: fluid state('L': liquid, 'V': vapor) 
% eos: equation of state (VR, VDW, RK, SRK, PR) 
% P,T: pressure(bar) and temperature(K) 
% Tc,Pc: critical T(K) and P(bar) 
% w- acentric factor 
% output: 
% phig: fugacity coefficient 
% f: fugacity (bar) 

% Tr and Pr (reduced T and P) 
Tr = T/Tc; Pr = P/Pc; R = 83.14; % cm^3*bar/mol/K 
eos = upper(eos); 
switch eos    
case 'VDW', al = 1; sm = 0; ep = 0; om = 0.125; ps = 0.42188; kappa = 0;    
case 'RK', al = 1./sqrt(Tr); sm = 1; ep = 0; om = 0.08664; ps = 0.42748;     
case 'SRK'         
kappa = 0.480 + 1.574*w - 0.176*w^2;         
al = (1 + kappa*(1-sqrt(Tr))).^2; sm = 1; ep = 0; om = 0.08664; ps = 0.42748;    
otherwise % PR or VR         
kappa = 0.37464 + 1.54226*w - 0.26992*w^2;         
al = (1 + kappa*(1-sqrt(Tr))).^2; sm = 1+sqrt(2); ep = 1-sqrt(2);  
om = 0.0778; ps = 0.45724;  
end 
% compressibility factor (Z) 
state = upper(state); beta = om*Pr./Tr; q = ps*al./(om*Tr); % beta and q 
if strcmp(eos,'VR') % virial EOS: vapor phase    
B0 = 0.083 - 0.422./(Tr.^1.6); B1 = 0.139 - 0.172./(Tr.^4.2);     
B = R*Tc.*(B0 + w.*B1)./Pc; Z = 1 + B.*P./(R*T); 
else    
fV = @(Z) 1+beta-q*beta.*(Z-beta)./((Z+ep*beta).*(Z+sm*beta)) - Z;     
fL = @(Z) beta+(Z+ep*beta).*(Z+sm*beta).*(1+beta-Z)./(q.*beta) - Z;     
switch state        
case 'V', Z = fzero(fV, 1);         
case 'L', Z = fzero(fL, beta);    
end 
end 
V = Z*R*T/P; % cm^3/mol 
% fugacity coefficient 
a = ps*al*R^2*Tc.^2/Pc; b = om*R*Tc./Pc; % a, b 
qi = a/(b*R*T); Bd = b*P/(R*T); % A, B 
if strcmp(eos,'VR'), phig = exp(Pr*(B0 + w*B1)./Tr); % virial EOS: vapor phase 
elseif strcmp(eos,'VDW'), phig = exp(Z - 1 - log(Z.*(1 - b/V)) - a./(R*T*V)); 
elseif strcmp(eos,'RK'), phig = exp(Z-1 - log(Z.*(1 - b/V)) - (a/(b*R*T)) * log(1+b/V)); 
else % SRK or PR    
phig = exp(Z-1 - log(Z-Bd) - (qi/(sm-ep))*log((Z + sm*Bd)/(Z + ep*Bd))); 
end 
f = phig*P; % fugacity (bar) 
end