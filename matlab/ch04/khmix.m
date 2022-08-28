function [Z H] = khmix(x,P,T,state,eos, Pc,Tc,w,k,Afi) 
% Estimation of enthalpy of mixture 
% Inputs: 
% x: mole fractions of all components (column vector) 
% P, T: pressure(Pa) and temperature(K) 
% state: fluid state('L': liquid, 'V': vapor) 
% eos: equation of state('RK', 'SRK', or 'PR') 
% Pc, Tc: critical P(Pa) and T(K) (column vector) 
% w- acentric factors of all components (column vector) 
% k: binary interaction parameter matrix(n x n symmetric) 
% Afi: coefficients of ideal gas heat capacity relation(J/mol/K) 
% Outputs: 
% Z: compressibility factor 
% H: enthalpy of mixture (J/mol) (1btu/lbmole = 2.326 J/mol) 
R = 8.314; % gas constant: m^3 Pa/(mol K) = J/mol-K 
Tr = T./Tc; Pr = P./Pc; nc = length(x); eos = upper(eos); state = upper(state); 
switch eos 
case{'RK'}        
ai = sqrt(0.4278./(Pc.*Tr.^2.5)); bi = 0.0867./(Pc.*Tr);        
A = sum(x.*ai); B = sum(x.*bi); Z = roots([1 -1 B*P*(A^2/B-B*P-1) -A^2*(B*P)^2/B]);    
case{'SRK'}        
mx = 0.48+1.574*w-0.176*w.^2; al = (1+mx.*(1-sqrt(Tr))).^2;        
ai = 0.42747*al.*Pr./(Tr.^2); bi = 0.08664*Pr./Tr;  
am = sqrt(ai'*ai).*(1-k);        
A = x*am*x'; B = sum(x.*bi); Z = roots([1 -1 A-B-B^2 -A*B]);    
case{'PR'}        
mx = 0.37464+1.54226*w-0.26992*w.^2; al = (1+mx.*(1-sqrt(Tr))).^2;        
ai = 0.45723553*al.*Pr./(Tr.^2); bi = 0.0777961*Pr./Tr;  
am = sqrt(ai'*ai).*(1-k);        
A = x*am*x'; B = sum(x.*bi); Z = roots([1 B-1 A-3*B^2-2*B B^3+B^2-A*B]); 
end 
iz = abs(imag(Z)); Z(and(iz>0,iz<=1e-6)) = real(Z(and(iz>0,iz<=1e-6))); 
for i = 1:length(Z), zind(i) = isreal(Z(i)); end 
Z = Z(zind); 
if state == 'L', Z = min(Z); else, Z = max(Z); end 
V = R*T*Z/P; % m^3/mol 
% Hv0 = int(Cp) 
Tf = (T-273.15)*1.8 + 32; % T: K->F 
Hv0 = x*(Afi(:,1)*Tf + Afi(:,2)*Tf^2/2 + Afi(:,3)*Tf^3/3 + Afi(:,4)*Tf^4/4 + Afi(:,5)*Tf^5/5); 
Hv0 = Hv0*2.326; % Btu/lbmole -> J/mol 
% compute enthalpy 
switch eos    
case{'RK'}, H = Hv0 + R*T*(Z - 1 - 3*(A^2)*log(1 + B*P/Z)/(2*B));    
case{'SRK'}        
hsum = 0;        
for i = 1:nc            
for j = 1:nc                
hsum = hsum + x(i)*x(j)*am(i,j)*(1 - mx(i)*sqrt(Tr(i))/(2*sqrt(al(i))) -...                    
mx(j)*sqrt(Tr(j))/(2*sqrt(al(j))));            
end        
end        
H = Hv0 + R*T*(Z - 1 - log((Z + B)/Z)*hsum/B);    
case{'PR'}        
hsum = 0;        
for i = 1:nc            
for j = 1:nc                
hsum = hsum + x(i)*x(j)*am(i,j)*(1 - mx(i)*sqrt(Tr(i))/(2*sqrt(al(i))) -...                    
mx(j)*sqrt(Tr(j))/(2*sqrt(al(j))));            
end        
end        
H = Hv0 + R*T*(Z - 1 - log((Z + B)/Z)*hsum/B); 
end 
end 