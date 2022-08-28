% rxnf.m
Cas = 0.2; R = 0.5; De = 0.1; k1a = 6.4; rspan = [0 R]; critN = 1e-6; errN = 1; Ca1 =1e-2; Ca2 = 3e-2;
while errN > critN
Cam = (Ca1+Ca2)/2;
[r x1] = ode45(@rxf,rspan,[0,0,Ca1],[],Cas,R,De,k1a);
[r x2] = ode45(@rxf,rspan,[0,0,Ca2],[],Cas,R,De,k1a);
[r xm] = ode45(@rxf,rspan,[0,0,Cam],[],Cas,R,De,k1a);
if (x1(end,3)-Cas)*(xm(end,3)-Cas) < 0, Ca2 = Cam;
else, Ca1 = Cam; end
errN = abs(Ca1 - Ca2);
end
Na = xm(:,1)./r.^2; effc = xm(end,2); Ca = xm(:,3);
ephi = R*sqrt(k1a/De); effa = 3*(ephi*coth(ephi)-1)/ephi^2;
fprintf('Effectiveness factor at r=R (calculated): %7.5f\n', effc);
fprintf('Effectiveness factor at r=R (analytic): %7.5f\n', effa);
plot(r,Ca), xlabel('r(cm)'), ylabel('C_A(gmol/cm^3)'), grid

function dxdr = rxf(r,x,Cas,R,De,k1a)
% x1=Nar^2, x2=eta, x3=Ca
if r == 0, Na = 0; else Na = x(1)/r^2; end
dxdr = [-k1a*x(3)*r^2; 3*x(3)*r^2/(Cas*R^3); -Na/De];
end
