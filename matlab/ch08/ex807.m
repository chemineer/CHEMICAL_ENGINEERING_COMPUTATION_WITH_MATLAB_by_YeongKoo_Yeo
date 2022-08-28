% bisecQL.m: calculation of heat loss from pipe flange
clear all;
T0 = 260; Ta = 60; R1 = 0.0833; R2 = 0.25; B = 0.5/12; k = 133; h = 3;
qa = 200; qb = 800; % assume initial range for heat loss flux
feps = 1e-3; ferr = 1e3; iter = 1;
while ferr > feps
qm = (qa+qb)/2; % midpoint of the assumed temperature range
[r,ya] = ode45(@funQL,[R1 R2],[T0 qa],[],Ta,k,h,B);
[r,yb] = ode45(@funQL,[R1 R2],[T0 qb],[],Ta,k,h,B);
[r,ym] = ode45(@funQL,[R1 R2],[T0 qm],[],Ta,k,h,B);
fa = ya(end,2) - R2*h*(ya(end,1)-Ta); fb = yb(end,2) - R2*h*(yb(end,1)-Ta);
fm = ym(end,2) - R2*h*(ym(end,1)-Ta);
if (fa*fm) < 0, qb = qm;
else, qa = qm; end
ferr = abs(qb - qa); iter = iter+1;
end
Qr = ym(:,2)./r; T = ym(:,1); qrate = 2*pi*ym(:,2)*B;
subplot(1,2,1), plot(r,qrate), xlabel('r(ft)'), ylabel('q(Btu/hr)'), grid,
axis tight;
subplot(1,2,2), plot(r,T), xlabel('r(ft)'), ylabel('T(F)'), grid;
fprintf('Number of iterations: %d, heat transfer rate at r=R2 (Btu/h): %g\n', iter, qrate(end))
fprintf('T at r=R2 (K): %g\n', T(end))
function dy = funQL(r,y,Ta,k,h,B)
q = y(2)/r; dy(1) = -q/k; dy(2) = -h*r*(y(1)-Ta)/B; dy = dy';
end
