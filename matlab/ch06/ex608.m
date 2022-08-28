% pHcstr.m: a neutralization cstr 
clear all; 
x1i = 1.2e-3; x2i = 2e-3; x3i = 2.5e-3; Kd = 1e-7; Ke = 1e-14; % data 
Fa = 0.01667; Fb = 2.333e-3; V = 2.5; av = Fa/V; bv = Fb/V; % data 
xi = [x1i x2i x3i];  
% define ODE system 
dx = @(t,x) [av*(xi(1) - x(1)) - bv*x(1); bv*(xi(2) - x(2)) - av*x(2); bv*(xi(3) - x(3)) - av*x(3)]; 
tspan = [0 400]; [t x] = ode45(dx,tspan,xi); % Solve ODEs using ode45 
x1 = x(:,1); x2 = x(:,2); x3 = x(:,3); 
% Determine [H+] 
Ken = Ke*ones(length(t),1); Kdn = Kd*ones(length(t),1);  
f = @(h) h+x2+x3-x1-Ken./h-x3./(1+Kdn.*h./Ken); 
H0 = 1e-3*ones(length(t),1);  % Initial guesses 
Hp = fsolve(f,H0); pH = -log10(Hp); 
% Plot concentrations and pH 
subplot(1,2,1), plot(t,x1,t,x2,':',t,x3,'--'), grid, xlabel('t(sec)'), 
ylabel('x_i(mol/liter)'), legend('x_1','x_2','x_3') 
subplot(1,2,2), plot(t,pH), grid, xlabel('t(sec)'), ylabel('pH') 