% denvir.m 
B = -0.0361; C = 2.7047e-3; D = -4.4944e-4; % virial constants 
R = 0.08206; T = 200; n = 300; P = linspace(1,30,n); x0 = 0.5*ones(1,n); 
f = @(x) 1 + B*x + C*x.^2 + D*x.^3 - P./(x*R*T); % define nonlinear equation 
rho = fsolve(f,x0); plot(P,rho), xlabel('P(atm)'), ylabel('Density(mol/liter)'), grid 