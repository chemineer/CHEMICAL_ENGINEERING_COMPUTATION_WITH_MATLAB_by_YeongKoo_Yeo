function [Z V] = virialEOS(P,T,Pc,Tc,w) 
% Estimation of compressibility factor and molar volume  
% at given T and P using the virial equation of state 
% inputs: 
% P,Pc: pressure and critical pressure (atm) 
% T,Tc: temperature and critical temperature (K) 
% w: acentric factor 
% outputs: 
% Z: compressibility factor 
% V: molar volume 
R = 0.08206; % atm-liter/(gmol-K) 
Tr = T./Tc; Pr = P./Pc; B0 = 0.083 - 0.422./(Tr.^1.6); B1 = 0.139 - 0.172./(Tr.^4.2); 
B = R*Tc.*(B0 + w.*B1)./Pc; Z = 1 + B.*P./(R*T); V = R*Z.*T./P; 
end 