function [Z V dH dS] = deptfun(state,eos,T,P,Tc,Pc,w) 
% Calculation of departure functions using the virial and cubic EOS 
% inputs 
% state: fluid state (liquid: L, vapor: V) 
% eos: type of the equation of state (VR, VDW, RK, SRK, PR) 
% T,P: temperature (K) and pressure (bar) 
% Tc,Pc: critical temperature (K) and critical pressure (bar) 
% w- acentric factor 
% outputs: 
% Z: compressibility factor 
% V: molar volume 
% dH: enthalpy departure 
% dS: entropy departure 
% Tr and Pr (reduced T and P) 
Tr = T/Tc; Pr = P/Pc; R = 83.14; % cm^3*bar/mol/K 
eos = upper(eos); 
switch eos    
case 'VDW', al = 1; sm = 0; ep = 0; om = 0.125; ps = 0.42188; kappa = 0;    
case 'RK', al = 1./sqrt(Tr); sm = 1; ep = 0; om = 0.08664; ps = 0.42748;     
case 'SRK', kappa = 0.480 + 1.574*w - 0.176*w^2;          
al = (1 + kappa*(1-sqrt(Tr))).^2; sm = 1; ep = 0; om = 0.08664; ps = 0.42748;    
otherwise % PR, or VR         
kappa = 0.37464 + 1.54226*w - 0.26992*w^2; al = (1 + kappa* (1-sqrt(Tr))).^2;         
sm = 1+sqrt(2); ep = 1-sqrt(2); om = 0.0778; ps = 0.45724;  
end 
% compressibility factor (Z) 
state = upper(state); beta = om*Pr./Tr; q = ps*al./(om*Tr); % beta and q 
if strcmp(eos,'VR') % virial (VR) EOS: vapor phase    
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
% departure function 
a = ps*al*R^2*Tc.^2/Pc; b = om*R*Tc./Pc; % a, b 
Ad = a*P/(R^2*T^2); Bd = b*P/(R*T); % A, B 
if strcmp(eos,'VR') % virial EOS: vapor phase    
dH = -Pr*(1.0972/Tr^2.6-0.083/Tr+w*(0.8944/Tr^5.2-0.139/Tr))*R*T;     
dS = -Pr*(0.675/Tr^2.6 + w*0.722/Tr^5.2)*R; 
elseif strcmp(eos,'RK') % Redlich-Kwong EOS    
dH = (Z - 1 - 1.5*ps/(om*Tr^1.5) * log(1 + b/V))*R*T;     
dS = (log(Z-Bd) - 0.5*ps/(om*Tr^1.5) * log(1 + b/V))*R; 
elseif strcmp(eos,'VDW') % van der Waals EOS    
dH = (Z - 1 - 3.375*Bd/Tr/Z)*R*T;     
dS = R*log(Z-Bd); 
elseif strcmp(eos,'SRK') % Soave-Redlich-Kwong EOS    
dH = (Z - 1 + (-kappa*sqrt(Tr)/(1 + kappa*(1-sqrt(Tr))) - 1)*(Ad/Bd) *log(1 + Bd/Z))*R*T;    
dS = (log(Z - Bd) - kappa*sqrt(Tr)/(1 + kappa*(1-sqrt(Tr)))*(Ad/Bd) *log(1 + Bd/Z))*R; 
else % Peng-Robinson EOS    
dH = (Z-1 -(Ad/(Bd*sqrt(8)))*(1+ kappa*sqrt(Tr)/sqrt(al)) *log((Z + sm*Bd)/(Z + ep*Bd)))*R*T;    
dS = (log(Z-Bd) - (Ad/(Bd*sqrt(8)))*(kappa*sqrt(Tr)/sqrt(al)) *log((Z + sm*Bd)/(Z + ep*Bd)))*R; 
end 
% R = 83.14 cm^3*bar/mol/K = 8.314 J/mol/K 
dH = dH/10; dS = dS/10; % dH:J/mol, dS:J/mol/K 
end 