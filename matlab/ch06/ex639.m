% micbrx.m
clear all;
% Define data structure (mbdat)
mbdat.Y = 0.4; mbdat.D = 0.202; mbdat.alp1 = 2.2; mbdat.alp2 = 0.2;
mbdat.Pf = 50; mbdat.mum = 0.48; mbdat.K1 = 0.04545; mbdat.Km = 1.2; mbdat.Sf = 20;
% Solve ODEs
x0 = 1; S0 = 50; P0 = 0;  % Initial values
tint = [0 120]; z0 = [x0 S0 P0]; [t z] = ode45(@mbrxn,tint,z0,[],mbdat);
x = z(:,1); S = z(:,2); P = z(:,3);
subplot(1,2,1), plot(t,x,t,P,'--'), grid, xlabel('t(h)'), ylabel('x and P (g/l)'), legend('x','P','location','best')
subplot(1,2,2), plot(t,S), grid, xlabel('t(h)'), ylabel('S (g/l)')
% Determine steady-state values
xs0 = 5; Ss0 = 5; Ps0 = 10; zs0 = [xs0 Ss0 Ps0]; zs = fsolve(@mbsrxn,zs0,[],mbdat);
fprintf('At steady-state, x = %g, S = %g, P = %g\n', zs(1), zs(2), zs(3));
function dz = mbrxn(t,z,mbdat)
% z(1)=x, z(2)=S, z(3)=P
Y = mbdat.Y; D = mbdat.D; alp1 = mbdat.alp1; alp2 = mbdat.alp2;
Pf = mbdat.Pf; mum = mbdat.mum; K1 = mbdat.K1; Km = mbdat.Km; Sf = mbdat.Sf;
mu = mum*z(2)*(1-z(3)/Pf)/(Km+z(2)+K1*(z(2))^2);
dz = [(mu - D)*z(1); D*(Sf - z(2)) - mu*z(1)/Y; -D*z(3) + (alp1*mu + alp2)*z(1)];
end
function f = mbsrxn(z,mbdat)
% z(1)=x, z(2)=S, z(3)=P
Y = mbdat.Y; D = mbdat.D; alp1 = mbdat.alp1; alp2 = mbdat.alp2;
Pf = mbdat.Pf; mum = mbdat.mum; K1 = mbdat.K1; Km = mbdat.Km; Sf = mbdat.Sf;
mu = mum*z(2)*(1-z(3)/Pf)/(Km+z(2)+K1*(z(2))^2);
f = [(mu - D)*z(1); D*(Sf - z(2)) - mu*z(1)/Y; -D*z(3) + (alp1*mu + alp2)*z(1)];
end
