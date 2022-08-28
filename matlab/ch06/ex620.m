% pfrnon.m: non-isothermal PFR using the method of line (MoL) 
clear all; 
% Data and parameters  
n = 50; pf.L = 2; pf.v = 0.4; pf.n = n; pf.Cf = 1; pf.Tf = 450; pf.E = 6e4;  
pf.dEr = -1e4; pf.dH = -1e5; pf.rCp = 800; pf.T1 = 450; pf.k1 = 0.2; pf.Kr1 = 1; 
pf.R = 8.314;  
% Initialization 
Cf = pf.Cf; Tf = pf.Tf; h = pf.L/n; w = [0:n]*h; w = w(2:end); 
Z0 = [ones(1,n)*Cf; ones(1,n)*Tf]; Z0 = reshape(Z0,2*n,1); tspan = [0:0.1:10]; 
% Employ stiff system solver  
[t Z] = ode15s(@dfr,tspan,Z0,[],pf); 
C = Z(:,1:2:end); T = Z(:,2:2:end); % Concentration and temperature 
% Display results 
subplot(1,2,1), mesh(w,t,C); xlabel('x(m)'), ylabel('t(min)'), zlabel('C(mol/m^3)'), colormap(gray) 
subplot(1,2,2), mesh(w,t,T); xlabel('x(m)'), ylabel('t(min)'), zlabel('T(K)'), colormap(gray) 

function dZ = dfr(t,Z,pf) 
% Non-isothermal PFR 
% Retrieve data 
L = pf.L; v = pf.v; n = pf.n; Cf = pf.Cf; Tf = pf.Tf; E = pf.E; dEr = pf.dEr; 
T1 = pf.T1; k1 = pf.k1; Kr1 = pf.Kr1; dH = pf.dH; R = pf.R; rCp = pf.rCp; 
% Initialization 
h = L/n; Z = reshape(Z,2,n); C = Z(1,:); T = Z(2,:); dC = zeros(n,1); dT = zeros(n,1);  
% Difference equations 
for i = 1:n    
k = k1*exp(-(E/R)*(1/T(i) - 1/T1)); Kr = Kr1*exp(-(dEr/R)*(1/T(i) - 1/T1));     
rx = k*C(i)/sqrt(1 + Kr*C(i)^2);    
if i == 1, s = (v/h)*(C(i) - Cf); d = (v/h)*(T(i) - Tf);    
else, s = (v/h)*(C(i) - C(i-1)); d = (v/h)*(T(i) - T(i-1)); end    
dC(i) = - s - rx; dT(i) = - d + (-dH)*rx/rCp; % difference model equation 
end 
dZ = reshape([dC'; dT'],2*n,1); 
end 