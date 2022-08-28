% csthcont.m
% z(1) = T, z(2) = T1, z(3) = Tt, z(4) = Ce
rho = 980; V = 2.5; Cp = 1.6; w = 250; taut = 3.6; theta = 1.2; Tr = 75;
Ti = 50; tauI = 1.8; K = 1/(Cp*w); tau = rho*V/w; Qs = (Tr-Ti)/K;
z0 = [Tr Tr Tr 0]; % Integral time interval and initial conditions
tv = [0 60]; Kc = 0;
[t z] = ode45(@htf,tv,z0,[],tau,taut,tauI,K,Kc,theta,Qs,Ti,Tr);
T = z(:,1); T1 = z(:,2); Tt = z(:,3);
plot(t,T,t,Tt,':',t,T1,'.-'), xlabel('t(min)'), ylabel('T(C)')
legend('Tank','Thermocouple','Outlet','location','best')

function dz = htf(t,z,tau,taut,tauI,K,Kc,theta,Qs,Ti,Tr)
% z(1) = T, z(2) = T1, z(3) = Tt, z(4) = Ce
if t < 10, Ti = 50; else Ti = 30; end
Q = Qs + (Kc/tauI)*z(4) + Kc*(Tr - z(3));
dz(1,1) = (Ti - z(1))/tau + K*Q/tau;
dz(2,1) = 2*(z(1) - z(2) - (theta/2)*dz(1,1))/theta;
dz(3,1) = (z(2) - z(3))/taut; dz(4,1) = Tr - z(3);
end
