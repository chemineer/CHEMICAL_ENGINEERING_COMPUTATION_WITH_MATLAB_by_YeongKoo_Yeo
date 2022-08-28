% Liquid Flow in a Pipeline
P1=205; P2=125; dP=P2-P1; dz=-125; rho=1000; g=9.81; Ws=1.2e6; 
mdot = -Ws/(dP/rho + g*dz) 