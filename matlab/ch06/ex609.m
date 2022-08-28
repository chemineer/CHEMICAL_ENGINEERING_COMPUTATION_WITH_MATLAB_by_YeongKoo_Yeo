% plcstr.m 
Ca0 = 0.4; k0 = 460; E = 1380; tau = 0.18; Tc = 298.15; kappa = 78; Cp = 32; % data
dy = @(t,x) [(Ca0 - x(1))/tau - (-(-k0*exp(-E/x(2))*x(1))); % define differential equations     
(-(-151080 + 2*(x(2) - 298.15))/Cp)*(-(-k0*exp(-E/x(2))*x(1))/Ca0) - (1+kappa)*(x(2) - Tc)/tau ]; 
tspan = [0 1]; y0 = [0.1 300]; [t y] = ode45(dy, tspan, y0); % solve ODE by ode45 
y(:,2) = y(:,2) - 273.15;  % K -> deg.C 
subplot(1,2,1), plot(t,y(:,1)), xlabel('t(hr)'), ylabel('C_A(mol/cm^3)') 
subplot(1,2,2), plot(t,y(:,2)), xlabel('t(hr)'), ylabel('T(deg.C)') 