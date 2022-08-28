% pfrLHmodel.m: PFR reaction calculation by Langmuir-Hinshelwood model
clear all; 
q  =  0.12; Ac  =  0.26; k  =  2.1; kr  =  0.98; Ca0  =  1; % data 
% (1) Profiles of concentration and conversion 
zspan  =  [0 0.5]; 
dC  =  @(z,Ca) -(Ac/q)*k*Ca/sqrt(1 + kr*Ca^2); 
[z Ca]  =  ode45(dC, zspan, Ca0); x  =  (Ca0 - Ca)/Ca0; 
subplot(1,2,1), plot(z,Ca), grid, xlabel('z(m)'), ylabel('C_A(mol/m^3)') 
subplot(1,2,2), plot(z,x), grid, xlabel('z(m)'), ylabel('x(conversion)') 
% (2) Calculation of reactor volume 
V0  =  0; xd  =  0.8; xspan  =  [0 xd]; % conversion of A  =  80% 
dV  =  @(x,V) (q/k)*sqrt(1 + kr*Ca0^2*(1-x)^2)/(1-x); 
[x V]  =  ode45(dV, xspan, V0); 
fprintf('Reactor volume for 80 percent conversion  =  %g m^3\n', V(end))