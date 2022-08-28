% bisecQT.m: calculation of heat flux and temperature in a wire
clear all;
T1 = 288.15; R1 = 0.004; k = 5; Ir = 400/(pi*R1^2);
T0a = 300; T0b = 400; % guess initial temperature
Teps = 1e-3; Terr = 1e3; iter = 1;
while Terr > Teps
T0m = (T0a+T0b)/2;
[r,za] = ode45(@funQT,[0 R1],[T0a 0],[],k,R1,Ir);
[r,zb] = ode45(@funQT,[0 R1],[T0b 0],[],k,R1,Ir);
[r,zm] = ode45(@funQT,[0 R1],[T0m 0],[],k,R1,Ir);
a = za(end,1); b = zb(end,1); m = zm(end,1);
if (m-T1) > 0, T0b = T0m;
else, T0a = T0m; end
Terr = abs(T0b - T0a); iter = iter+1;
end
Qr = zm(:,2)./r; T = zm(:,1);
subplot(1,2,1), plot(r,Qr), xlabel('r(m)'), ylabel('Qr(W/m^2)'), grid;
subplot(1,2,2), plot(r,T), xlabel('r(m)'), ylabel('T(K)'), grid;
fprintf('Number of iterations: %d, Qr at r=R1 (W/m^2): %g, ', iter, Qr(end))
fprintf('T at r=0 (K): %g\n', T(1))

function dz = funQT(r,z,k,R1,Ir) 
if r > 0, Qr = z(2)/r; 
else, Qr = 0; end 
dz(1) = -Qr/k; dz(2) = Ir^2*r/(1.4e5*exp(0.0035*z(1))); dz = dz'; 
end 