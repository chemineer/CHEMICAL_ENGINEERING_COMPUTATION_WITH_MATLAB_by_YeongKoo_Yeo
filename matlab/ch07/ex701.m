% xaz.m
Dab = 1.991e-5; T = 328.5; R = 8.314; P = 99.4; Pa0 = 68.4; Tf = 295; z2 = 0.238; z1 = 0; zv = [z1 z2];
C = P/(R*T); x0 = Pa0/P; critN = 1e-11; errN = 1; Na1 = 3.5e-6; Na2 = 3.6e-6;
while errN > critN
Nam = (Na1+Na2)/2;
[z x1] = ode45(@dxz,zv,x0,[],Na1,Dab,C);
[z x2] = ode45(@dxz,zv,x0,[],Na2,Dab,C);
[z xm] = ode45(@dxz,zv,x0,[],Nam,Dab,C);
if x1(end)*xm(end) < 0, Na2 = Nam;
else, Na1 = Nam; end
errN = abs(Na1 - Na2);
end
xblm = x0/(log(1/(1-x0))); Nanal = Dab*C*x0/((z2-z1)*xblm);
fprintf('Estimated Nab = %e, xA = %8.6f\n',Nam, xm(end));
fprintf('Analytic Nab = %e\n',Nanal);
plot(z,xm), xlabel('z(m)'), ylabel('x_A'), grid, axis tight

function dx = dxz(z,x,Na,Dab,C)
dx = -(1-x)*Na/Dab/C;
end
