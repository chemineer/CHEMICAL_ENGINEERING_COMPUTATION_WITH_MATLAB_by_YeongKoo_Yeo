% so3rxn.m: oxidation of SO2 using fixed-bed reactor 
clear all; 
% Data 
sa.Ta = 700; sa.T0 = 750; sa.Pt0 = 202650; % initial temperature and pressure 
sa.rho0 = 0.866; sa.rhob = 542; % density 
sa.ya0 = 0.1; sa.yb0 = 0.11; sa.yc0 = 0.79; % a:SO2, b:O2, c:N2, d:SO3  
sa.Ft0 = 0.02153; sa.G = 0.433; % feed rate(Ft0) and mass velocity(G) 
sa.epn = -0.05; sa.phi = 0.45; sa.mu = 3.72e-5; % parameters and viscosity 
sa.D = 0.0453; sa.Dp = 4.57e-3; % tube diameter(D) and particle diameter(Dp) 
sa.U = 17; % overall heat transfer coefficient(J/s/m^2/K) 
% Solve differential equations: x = z(1); T = z(2); P = z(3); 
wspan = [0 4]; z0 = [0, sa.T0, sa.Pt0]; [w z] = ode15s(@so3de,wspan,z0,[],sa); 
% Plot results 
x = z(:,1); T = z(:,2); P = z(:,3); 
subplot(2,2,1), plot(w,x), grid, xlabel('W(kg)'), ylabel('X') 
subplot(2,2,2), plot(w,T), grid, xlabel('W(kg)'), ylabel('T(K)') 
subplot(2,2,3), plot(w,P), grid, xlabel('W(kg)'), ylabel('P(Pa)') 