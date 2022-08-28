% chemostat.m: calculation of a chemostat (continuous stirred biological reactor)
clear all;
D = 0.1; Sf = 5; y1 = 0.8; y2 = 0.7; mum = 0.6; Km = 0.28; % data
X0 = 0.03; S0 = 5; P0 = 0; % initial conditions (z(1) = x, z(2) = S, z(3) = P)
dz = @(t,z) [-D*z(1) + y1*mum*z(2)*z(1)/(Km + z(2));
- mum*z(2)*z(1)/(Km + z(2)) + D*(Sf - z(2));
y2*mum*z(2)*z(1)/(Km + z(2)) - D*z(3)];
tint = [0 30]; z0 = [X0, S0, P0];
[t z] = ode15s(dz,tint,z0); % z(1) = x, z(2) = S, z(3) = P
x = z(:,1); S = z(:,2); P = z(:,3);
plot(t,x,t,S,'--',t,P,':'), xlabel('t(hr)')
ylabel('Concentration(g/ml)'), legend('x','S','P','location','best')
