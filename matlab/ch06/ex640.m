% predprey.m: predator-prey model (E.Coli - Amoeba)
clear all;
pdat.Fv = 0.0625; pdat.S0 = 0.5; pdat.m1 = 0.25; pdat.m2 = 0.24; % data
pdat.K1 = 5e-4; pdat.K2 = 4e8; pdat.c1 = 3.3e-10; pdat.c2 = 1.4e3;
N10 = 1.3e9; N20 = 4e5; % initial conditions
x0 = [pdat.S0, N10, N20]; tspan = [0 1000];
[t x] = ode15s(@pred,tspan,x0,[],pdat);
plot(t,log10(x(:,2)),t,log10(x(:,3)),'--'), xlabel('t(hr)'), ylabel('log_{10}(N_1) and log_{10}(N_2)')
legend('Bacteria(N_1)','Amoeba(N_2)','location','best')
function dx = pred(t,x,pdat)
S = x(1); N1 = x(2); N2 = x(3);
% Retrieve data
S0 = pdat.S0; Fv = pdat.Fv; c1 = pdat.c1; c2 = pdat.c2;
m1 = pdat.m1; m2 = pdat.m2; K1 = pdat.K1; K2 = pdat.K2;
% Define differential equations
dx = [Fv*(S0 - S) - (c1*m1*N1*S)/(K1 + S);
-Fv*N1 + (m1*N1*S)/(K1 + S) - (c2*m2*N1*N2)/(K2 + N1);
-Fv*N2 + (m2*N1*N2)/(K2 + N1)];
end
