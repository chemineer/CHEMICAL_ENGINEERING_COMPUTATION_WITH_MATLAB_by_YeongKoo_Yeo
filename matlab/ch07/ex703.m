% sphdf.m
Dab = 7.39e-6; T = 298.15; P = 1.01325e5; p0 = 133.32; r0 = 0.003;
xspan = [0 20]; critN = 1e-16; errN = 1; Na1 = 1e-12; Na2 = 2e-12; % initial guesses
while errN > critN % bisection method
Nam = (Na1+Na2)/2;
[x z1] = ode45(@spf,xspan,[Na1,p0],[],T,P,Dab,r0);
[x z2] = ode45(@spf,xspan,[Na2,p0],[],T,P,Dab,r0);
[x zm] = ode45(@spf,xspan,[Nam,p0],[],T,P,Dab,r0);
if z1(end,2)*zm(end,2) < 0, Na2 = Nam;
else, Na1 = Nam; end
errN = abs(Na1 - Na2);
end
r = r0 + x; Na = zm(:,1)./r.^2; pf = zm(end,2);
fprintf('Flux at r=r0: %7.5e, partial pressure at r=inf: %7.5f\n', Na(1), pf);
plot(r(1:25),Na(1:25)), xlabel('r(m)'), ylabel('N_A')

function dzdx = spf(x,z,T,P,Dab,r0)
R = 8314.34;
dzdx = [0; -R*T*z(1)*(1-z(2)/P)/(Dab*(x+r0)^2)];
end
