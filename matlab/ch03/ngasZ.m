function nz = ngasZ(T,P,Sg) 
% Estimates the compressibility factor Z of natural gases 
% input 
% T: temperature (F) (scalar or vector) 
% P: pressure (psia) 
% Sg: specific gravity of the natural gas 
% output 
% nz: estimated compressibility factor Z of natural gases 
P = P/1000; T = T + 460; 
A1 = 0.001946; A2 = -0.027635; A3 = 0.136315; A4 = -0.23849; A5 = 0.105168; A6 = 3.44e8;  
F1 = P.*(0.251*Sg-0.15) - 0.202*Sg + 1.106; den = 1 + A6*P.*10.^(1.785*Sg)./ (T.^3.825); 
F2 = 1.4*exp(-0.0054*(T-460)); F3 = A1*P.^5 + A2*P.^4 + A3*P.^3 + A4*P.^2 + A5*P; 
F4 = (0.154-0.152*Sg).*P.^(3.18*Sg-1).*exp(-0.5*P) - 0.02;  
F5 = 0.35*(0.6-Sg).*exp(-1.039*(P-1.8).^2); 
nz = F1.*(1./den + F2.*F3) + F4 + F5; 
end 