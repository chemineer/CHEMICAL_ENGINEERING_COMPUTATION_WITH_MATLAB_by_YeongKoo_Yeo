% absf.m
% x(1)=Ca2, x(2)=Ca3, ..., x(9)=Ca10
Dab = 1.5e-9; delt = 3e-4; delx = delt/10; vm = 0.6; cas = 0.03; kp = 0;
c0 = zeros(1,9); % initial concentration profile
zf = 1; zspan = [0 zf]; [z x] = ode45(@bsf,zspan,c0,[],Dab,kp,delx,delt,vm,cas);
nz = length(z); xc1 = cas*ones(1,nz); % Ca1
for k = 1:nz
if (4*x(k,9)<x(k,8)), xc11(k) = 0;
else, xc11(k) = (4*x(k,9) - x(k,8))/3; end
end
xc = [];
for k = 1:9, xc = [xc x(:,k)]; end
c = [xc1' xc xc11']; plot(z,c(:,3),z,c(:,5),':',z,c(:,7),'.-',z,c(:,9),'--')
xlabel('z(m)'), ylabel('C(kgmol/m^3)'), legend('C_{A3}','C_{A5}','C_{A7}','C_{A9}','location','best')

function dxdz = bsf(z,x,Dab,kp,delx,delt,vm,cas)
% x(1)=Ca2, x(2)=Ca3, ..., x(9)=Ca10
ca1 = cas;
if (4*x(9)<x(8)), ca11 = 0;
else, ca11 = (4*x(9) - x(8))/3; end
dxdz = [(Dab*(x(2)-2*x(1)+ca1)/(delx^2) - kp*x(1))/(vm*(1 - (1*delx/delt)^2));
(Dab*(x(3)-2*x(2)+x(1))/(delx^2) - kp*x(2))/(vm*(1 - (2*delx/delt)^2));
(Dab*(x(4)-2*x(3)+x(2))/(delx^2) - kp*x(3))/(vm*(1 - (3*delx/delt)^2));
(Dab*(x(5)-2*x(4)+x(3))/(delx^2) - kp*x(4))/(vm*(1 - (4*delx/delt)^2));
(Dab*(x(6)-2*x(5)+x(4))/(delx^2) - kp*x(5))/(vm*(1 - (5*delx/delt)^2));
(Dab*(x(7)-2*x(6)+x(5))/(delx^2) - kp*x(6))/(vm*(1 - (6*delx/delt)^2));
(Dab*(x(8)-2*x(7)+x(6))/(delx^2) - kp*x(7))/(vm*(1 - (7*delx/delt)^2));
(Dab*(x(9)-2*x(8)+x(7))/(delx^2) - kp*x(8))/(vm*(1 - (8*delx/delt)^2));
(Dab*(ca11-2*x(9)+x(8))/(delx^2) - kp*x(9))/(vm*(1 - (9*delx/delt)^2))];
end
