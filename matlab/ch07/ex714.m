% abscolm.m: absorption column film mass transfer
clear all;
x = linspace(0,0.01,100); % 100 values from 0 to 0.01
t = linspace(0,1,100); % 100 values from 0 to 1
m = 0;
res = pdepe(m,@abeqn,@abitcon,@abbncon,x,t);
u = res(:,:,1);
subplot(1,2,1), surf(x,t,u), colormap(gray), xlabel('\delta(m)'), ylabel('L(m)'), zlabel('C_a(kgmol/m^3)')
shading interp  % surface plot
subplot(1,2,2)
% plot Ca as a function of delta(m)
for k = 1:length(t), plot(x,u(k,:),'k'), xlabel('\delta(m)'), ylabel('C_a(kgmol/m^3)'), hold on, end
hold off

function [c,f,s] = abeqn(x,t,u,DuDx)
delt = 0.01; mu = 2.1e-5; c = 2*(1/2)*(x/delt - (1/2)*(x/delt)^2)/mu; f = DuDx; s = 0;
end
function u0 = abitcon(x)
u0 = 0;
end

function [pl,ql,pr,qr] = abbncon(xl,ul,xr,ur,t)
C0 = 0.1;
pl = 0; ql = 1; % y = 0
pr = ur - C0; qr = 0; % y = delta
end
