% htcylinder.m: heat transfer in cylindrical laminar flow
clear all;
R = 0.03; tmax = 0.5; m = 1;
x = linspace(0,R,40); % subdivision from 0 to R
t = linspace(0,tmax,20); % subdivision from 0 to tmax
res = pdepe(m,@hteqn,@htitcon,@htbncon,x,t);
u = res(:,:,1);
subplot(1,2,1), mesh(x,t,u), colormap(jet), xlabel('r(m)'), ylabel('L(m)'),
zlabel('T(K)')
shading interp  % surface plot
subplot(1,2,2)
% plot T as a function of r(m)
for k = 1:length(t), plot(x,u(k,:),'k'), xlabel('r(m)'), ylabel('TK)'), hold on, end
grid, hold off
function [c,f,s] = hteqn(x,t,u,DuDx)
Cp = 4.2; rho = 1e3; vk = 5; R = 0.03; c = rho*Cp*vk*(1 - (x/R)^2); f = DuDx; s = 0;
end
function u0 = htitcon(x)
T0 = 500; u0 = T0;
end

function [pl,ql,pr,qr] = htbncon(xl,ul,xr,ur,t)
pl = 0; ql = 1; % r = 0
k = 0.1; Q = 100; pr = -Q; qr = -k; % r = R
end
