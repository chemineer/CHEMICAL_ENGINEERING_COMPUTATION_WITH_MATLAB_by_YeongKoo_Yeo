function rhoL = gunyam(Tc,Pc,w,Mw,T) 
% saturated liquid molar volume by the Gunn and Yamada method 
% Tc: deg,C, Pc: bar, Mw: g/mol 
R = 0.08314; % liter*atm/(mol*K) 
Tr = (T+273.15)/(Tc+273.15);          
if Tr <= 0.8                   
Vr0 = 0.33593 - 0.33953*Tr + 1.51941*Tr^2 - 2.02512*Tr^3 + 1.11422*Tr^4;    
else 
Vr0 = 1 + 1.3*sqrt(1-Tr)*log10(1-Tr) - 0.50879*(1-Tr) - 0.91534*(1-Tr)^2;    
end 
Gam = 0.29607 - 0.09045*Tr - 0.04842*Tr^2; Vsc = (R*(Tc+273.15)/Pc)*(0.2920 - 0.0967*w); 
V = Vsc*Vr0.*(1 - w*Gam); rhoL = Mw./V; % density: kg/m^3  
end 