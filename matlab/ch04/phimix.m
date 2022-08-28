function [Z,V,phi] = phimix(ni,P,T,Pc,Tc,w,k,state,eos) 
% Estimates fugacity coefficients of all components in a mixture using 
% the cubic equation of state 
% Inputs: 
% ni: number of moles (or mole fractions) of each component (vector) 
% P, T: pressure(Pa) and temperature(K) 
% Pc, Tc: critical P(Pa) and T(K) of all components (vector) 
% w- acentric factors of all components (vector) 
% k: symmetric matrix of binary interaction parameters (n x n) 
% state: fluid state('L': liquid, 'V': vapor) 
% eos: equation of state('RK', 'SRK', or 'PR') 
% Outputs: 
% V: molar volume (m3/mol) 
% Z: compressibility factor 
% phi: fugacity coefficient vector 

ni = ni(:); Pc = Pc(:); Tc = Tc(:); w = w(:); % column vector 
x = ni/sum(ni); % mole fraction 
R = 8.314; % gas constant: m^3 Pa/(mol K) = J/mol-K 
Tref = 298.15; % reference T(K) 
Pref = 1e5; % reference P(Pa) 
Tr = T./Tc; eos = upper(eos); state = upper(state); 
switch eos    
case{'RK'}        
ep = 0; sm = 1; om = 0.08664; ps = 0.42748; al = 1./sqrt(Tr);    
case{'SRK'}        
ep = 0; sm = 1; om = 0.08664; ps = 0.42748;        
al = (1+(0.48+1.574*w-0.176*w.^2).*(1-sqrt(Tr))).^2;    
case{'PR'}        
ep = 1 - sqrt(2); sm = 1 + sqrt(2); om = 0.07780; ps = 0.45724; 
al = (1+(0.37464+1.54226*w-0.26992*w.^2).*(1-sqrt(Tr))).^2; 

end 
ai = ps*(R^2).*al.*(Tc.^2)./Pc; am = sqrt(ai*ai').*(1 - k); % nxn matrix 
a = x'*am*x; bi = om*R*Tc./Pc; b = x'*bi; beta = b*P/R/T; q = a/(b*R*T); 
% compressibility factor(Z) and molar volume (V) 
state = upper(state); 
c(1) = 1; c(2) = (sm +ep)*beta - (1+beta); 
c(3) = beta*(q + ep*sm*beta -(1+beta)*(sm+ep));  
c(4) = -beta^2*(q +(1+beta)*ep*sm); 
% Roots 
Z = roots(c); iz = abs(imag(Z)); Z(and(iz>0,iz<=1e-6)) =  real(Z(and(iz>0,iz<=1e-6))); 
for i = 1:length(Z), zind(i) = isreal(Z(i)); end 
Z = Z(zind); 
if state == 'L', Z = min(Z); else, Z = max(Z); end 
V = R*T*Z/P; 
% fugacity coefficients 
bara = (2*ni'*am - a*ones(1,length(ni)))'; barb = bi; 
phi = exp((Z - 1)*barb/b - log((V - b)*Z/V) + (a/(b*R*T))/(ep - sm)*...    
log((V + sm*b)/(V + ep*b))*(1 + bara/a - barb/b)); 
end 