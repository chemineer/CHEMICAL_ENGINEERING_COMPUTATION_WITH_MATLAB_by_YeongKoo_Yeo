% batPIcon.m: PI control for batch reactor
% Data
clear all;
A1 = 1.1; A2 = 172.2; E1 = 2.09e4; E2 = 4.18e4; R = 8.314;
dH1 = 4.18e4; dH2 = 8.36e4; rho = 1000; Cp = 1; Tc = 25;
Tsmax = 150; Tsmin = 70; Uj = 1.16; Ucmax = 4.42; Ucmin = 1.39;
AcV = 17; AjV = 30; Kc = 0.1; tauI = 360; us = 1; Tmax = 125;
dt = 0.1; tmax = 4000;
% Parameters
gam1 = dH1/(rho*Cp); gam2 = dH2/(rho*Cp);
a1 = (Uj*Tsmin*AjV + Ucmax*Tc*AcV)/(rho*Cp); a2 = -(Uj*AjV + Ucmax*AcV)/ (rho*Cp);
b1 = (Uj*AjV*(Tsmax - Tsmin) - (Ucmax - Ucmin)*AcV*Tc)/(rho*Cp);
b2 = (Ucmax - Ucmin)*AcV/(rho*Cp);
% Initialization
t = 0:dt:tmax+dt; n = length(t); T(1) = 25; Ca(1) = 1; Cb(1) = 0;
Td = 54 + 71*exp(-0.0025*t); % desired temperature trajectory
%Td(n+1) = 54 + 71*exp(-0.0025*(t(end)+dt));
er(1) = 100; erc(1) = 0; us = 1; u(1) = us;
Ts(1) = (Tsmax - Tsmin)*u(1) + Tsmin; Uc(1) = (Ucmin - Ucmax)*u(1) + Ucmax;
Fc(1) = (4550*(1/Uc(1) - 1/10.8))^(-1.25);
for k = 1:n-1
% Begin 4th-order RK method
T0 = T(k); Ca0 = Ca(k); Cb0 = Cb(k);
k1 = gam1*A1*exp(-E1/R/(273.15+T0))*Ca0^2 +...
gam2*A2*exp(-E2/R/(273.15+T0))*Cb0 + (a1+a2*T0) + (b1+b2*T0)*u(k);
k11 = -A1*exp(-E1/R/(273.15+T0))*Ca0^2;
k12 = A1*exp(-E1/R/(273.15+T0))*Ca0^2 - A2*exp(-E2/R/(273.15+T0))*Cb0;
T1 = T(k) + k1*dt/2; Ca1 = Ca(k) + k11*dt/2; Cb1 = Cb(k) + k12*dt/2;
k2 = gam1*A1*exp(-E1/R/(273.15+T1))*Ca1^2 +...
gam2*A2*exp(-E2/R/(273.15+T1))*Cb1 + (a1+a2*T1) + (b1+b2*T1)*u(k);
k21 = -A1*exp(-E1/R/(273.15+T1))*Ca1^2;
k22 = A1*exp(-E1/R/(273.15+T1))*Ca1^2 - A2*exp(-E2/R/(273.15+T1))*Cb1;
T2 = T(k) + k2*dt/2; Ca2 = Ca(k) + k21*dt/2; Cb2 = Cb(k) + k22*dt/2;
k3 = gam1*A1*exp(-E1/R/(273.15+T2))*Ca2^2 +...
gam2*A2*exp(-E2/R/(273.15+T2))*Cb2 + (a1+a2*T2) + (b1+b2*T2)*u(k);
k31 = -A1*exp(-E1/R/(273.15+T2))*Ca2^2;
k32 = A1*exp(-E1/R/(273.15+T2))*Ca2^2 - A2*exp(-E2/R/(273.15+T2))*Cb2;
T3 = T(k) + k3*dt/2; Ca3 = Ca(k) + k31*dt/2; Cb3 = Cb(k) + k32*dt/2;
k4 = gam1*A1*exp(-E1/R/(273.15+T3))*Ca3^2 +...
gam2*A2*exp(-E2/R/(273.15+T3))*Cb3 + (a1+a2*T3) + (b1+b2*T3)*u(k);
k41 = -A1*exp(-E1/R/(273.15+T3))*Ca3^2;
k42 = A1*exp(-E1/R/(273.15+T3))*Ca3^2 - A2*exp(-E2/R/(273.15+T3))*Cb3;
T(k+1) = T(k) + dt*(k1/6 + k2/3 + k3/3 + k4/6);
Ca(k+1) = Ca(k) + dt*(k11/6 + k21/3 + k31/3 + k41/6);
Cb(k+1) = Cb(k) + dt*(k12/6 + k22/3 + k32/3 + k42/6);
% End of 4th-order RK method
% PI control
erc(k+1) = erc(k) + er(k)*dt; er(k+1) = Td(k+1) - T(k+1);
u(k+1) = us + Kc*(er(k+1) + erc(k+1)/tauI);
if u(k+1) >= 1, u(k+1) = 1; elseif u(k+1) <= 0, u(k+1) = 0; end
% Calculate Ts, Uc and Fc
Ts(k+1) = (Tsmax - Tsmin)*u(k+1) + Tsmin;
Uc(k+1) = (Ucmin - Ucmax)*u(k+1) + Ucmax;
Fc(k+1) = (4550*(1/Uc(k+1) - 1/10.8))^(-1.25);
end
subplot(1,2,1), plot(t,Ca,'--',t,Cb), xlabel('t(sec)'), ylabel('C(kmol/m^3)')
legend('C_A(t)','C_B(t)'), axis([0 4000 0 1])
subplot(1,2,2), plot(t,Td,'--',t,T,t,Ts,'.-'), xlabel('t(sec)'), ylabel('Temp(deg.C)')
legend('Desired temp.','Reactor temp.','Steam temp.'), axis([0 4000 20 160])
