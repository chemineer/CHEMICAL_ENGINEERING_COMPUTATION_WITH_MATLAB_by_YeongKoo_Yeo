% httube.m: heat transfer in fluid flowing through a pipe
clear all;
% Data and parameters
L = 2; R = 0.05; % pipe length (L) and radius (R)
n = 20; h = R/n; r = [0:n]'*h; alpa = 1e-4; Ti = 300; Tb = 400; vmax = 0.5;
v = vmax*(1-(r/R).^2); % parabolic velocity profile
pf.r = r; pf.v = v; pf.n = n; pf.h = h; pf.alpa = alpa; pf.Tb = Tb;
% Solve PDE
T0 = ones(n,1)*Ti; opn = odeset('relTol',1e-12,'absTol',1e-10);
[x,T] = ode15s(@pflowht, [0,L],T0,opn,pf);
% Display results
na = length(x); Tw = Tb*ones(na,1); T = [T,Tw];
subplot(1,2,1), mesh(r,x,T), xlabel('r(m)'), ylabel('x(m)'),
zlabel('T(K)'), colormap(gray)
subplot(1,2,2), contourf(x,r,T'), xlabel('x(m)'), ylabel('r(m)'),
colormap(gray)
