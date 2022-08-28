% compdist.m
clear all;
% Define data structure (cmdat)
cmdat.F = 100; cmdat.z = 0.5; cmdat.R = 128.01; cmdat.V = 178.01; cmdat.mr = 10;
cmdat.mf = 10;
cmdat.ms = 10; cmdat.md = 100; cmdat.mb = 100; cmdat.alpa = 2;
% Solve ODEs: x(1)=xb, x(2)=xs, x(3)=xf, x(4)=xr, x(5)=xd
xb0 = 0; xs0 = 0; xf0 = cmdat.z; xr0 = 0; xd0 = 0; % Initial values
tint = [0 12]; x0 = [xb0 xs0 xf0 xr0 xd0];
[t x] = ode45(@cmeqn,tint,x0,[],cmdat); % solve differential eqn system
xb = x(:,1); xs = x(:,2); xf = x(:,3); xr = x(:,4); xd = x(:,5);
plot(t,xb,t,xs,'--',t,xf,':',t,xr,'.-',t,xd,'.'), grid, xlabel('t(min)')
ylabel('Composition'), legend('x_B','x_S','x_f','x_R','x_D','location','best')

function dx = cmeqn(t,x,cmdat) 
% x(1)=xb, x(2)=xs, x(3)=xf, x(4)=xr, x(5)=xd 
F = cmdat.F; z = cmdat.z; R = cmdat.R; V = cmdat.V; alpa = cmdat.alpa; 
mr = cmdat.mr; mf = cmdat.mf; ms = cmdat.ms; md = cmdat.md; mb = cmdat.mb; 
yb = alpa*x(1)/(1 + (alpa-1)*x(1)); ys = alpa*x(2)/(1 + (alpa-1)*x(2)); 
yf = alpa*x(3)/(1 + (alpa-1)*x(3)); yr = alpa*x(4)/(1 + (alpa-1)*x(4)); 
Lr = R; Ls = R + F; B = Ls - V; 
dx = [(Ls*x(2) - B*x(1) - V*yb)/mb;      
(Ls*(x(3) - x(2)) + V*(yb - ys))/ms;      
(Lr*(x(4) - x(3)) + F*(z - x(3)) + V*(ys - yf))/mf;      
(Lr*(x(5) - x(4)) + V*(yf - yr))/mr;      
V*(yr - x(5))/md]; 
end 
