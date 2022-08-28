% csthcontp.m
% z(1) = T, z(2) = T1, z(3) = Tt
rho = 980; V = 2.5; Cp = 1.6; w = 250; taut = 3.6; theta = 1.2; Tr = 75;
Ti = 50; Kc = 60; K = 1/(Cp*w); tau = rho*V/w; Qs = (Tr-Ti)/K;
tv = [0 60]; z0 = [Tr Tr Tr]; [t z] = ode45(@hpf,tv,z0,[],tau,taut,K,Kc,theta,Qs,Ti,Tr);
T = z(:,1); T1 = z(:,2); Tt = z(:,3);
plot(t,T,t,Tt,':',t,T1,'.-'), xlabel('t(min)'), ylabel('T(C)'), legend('Tank','Thermocouple','Outlet')

function dz = hpf (t,z,tau,taut,K,Kc,theta,Qs,Ti,Tr)
% z(1) = T, z(2) = T1, z(3) = Tt
if t < 10, Ti = 50; else Ti = 30; end
Q = Qs + Kc*(Tr - z(3));
dz(1,1) = (Ti - z(1))/tau + K*Q/tau;
dz(2,1) = 2*(z(1) - z(2) - (theta/2)*dz(1,1))/theta;
dz(3,1) = (z(2) - z(3))/taut;
end
