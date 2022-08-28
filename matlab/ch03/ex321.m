% Vapor Pressure of n-Propylbenzene
Tb = 159.22+273.15; nu = [2 1 1 5]; Dp = 0; GI = []; dB = [0.07545 -0.00227 0.11192 0.01653];  
T = 100 + 273.15; Pv = vpRM(T,Tb,nu,dB,GI,Dp); 

T = 200 + 273.15; Pv = vpRM(T,Tb,nu,dB,GI,Dp);