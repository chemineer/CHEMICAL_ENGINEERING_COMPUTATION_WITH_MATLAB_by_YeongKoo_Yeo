function Z = zLK(T,Tc,P,Pc,w) 
% zLK.m: calculates compressibility factor using the Lee-Kesler equation 
% Input: 
% T,Tc: temperature(K)
% P,Pc: pressure (MPa) 
% w: acentric factor 
% constants and parameters 
% Output:  
% Z: compressibility factor 

b = [0.1181193 0.2657280 0.1547900 0.0303230; 0.2026579 0.3315110 0.0276550 0.2034880]; 
c = [0.0236744 0.0186984 0.0000000 0.0427240; 0.0313385 0.0503618 0.0169010 0.0415770]; 
d = 1e-4*[0.155488 0.623689; 0.487360 0.0740336]; 
beta = [0.653920 1.22600]; gam = [0.060167 0.03754]; wr = 0.3978; 
Pr = P/Pc; Tr = T/Tc; % reduced pressure and temperature 
% Calculation of Z0 of simple fluid (ind = 1) 
ind = 1; [B C D] = compBCD(Tr,b,c,d,ind);  
Vr0 = Tr./Pr; Vr = fzero(@BWReq,Vr0,[],B,C,D,c,Tr,Pr,ind,beta,gam); Z0 = Pr*Vr./Tr; 
% Calculation of Zr of reference fluid (ind = 2) 
ind = 2; [B C D] = compBCD(Tr,b,c,d,ind); 
Vr = fzero(@BWReq,Vr,[],B,C,D,c,Tr,Pr,ind,beta,gam); Zr = Pr*Vr/Tr; 
Z1 = 1/wr * (Zr - Z0); 
Z = Z0 + w*Z1; % compressibility factor of real fluid 
fprintf('Compressibility factor Z = %f\n', Z); 
end  

function [B C D] = compBCD(Tr,b,c,d,ind) 
% ind 1: Simple fluids, ind 2: Reference fluids 
B = b(ind,1) - b(ind,2)/Tr - b(ind,3)/Tr^2 - b(ind,4)/Tr^3; 
C = c(ind,1) - c(ind,2)/Tr + c(ind,3)/Tr^3; D = d(ind,1) + d(ind,2)/Tr; 
end  

function f = BWReq(Vr,B,C,D,c,Tr,Pr,ind,beta,gam) 
f = 1 + B/Vr + C/Vr^2 + D/Vr^5 + c(ind, 4)/(Tr^3*Vr^2)*(beta(ind)+... 
gam(ind)/Vr^2)*exp(-gam(ind)/Vr^2) - Pr*Vr/Tr; 
end 