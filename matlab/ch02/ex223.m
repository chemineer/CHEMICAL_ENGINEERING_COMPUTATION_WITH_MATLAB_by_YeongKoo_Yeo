% N2denfit.m: density of N2 
clear all; T  =  200; R  =  0.08206; % gas constant 
P  =  [3.2 6.0 9.0 12.0 14.0 17.0 19.0 21.0]; % P(atm) 
rho  =  [0.1995 0.3825 0.5895 0.8108 0.9685 1.2257 1.4159 1.6281]; % rho(mol/l) 
Z  =  P./(R*rho*T); Vc  =  polyfit(rho, Z, 3) 
x  =  [0.15:0.01:1.8]; Zcal  =  polyval(Vc, x); Pcal  =  R*Zcal.*x*T; 
plot(x,Pcal,rho,P, 'o'), grid, xlabel ('\rho(mol/liter)'), ylabel ('P(atm)') 
legend('Fitting curve', 'Data','location','best') 