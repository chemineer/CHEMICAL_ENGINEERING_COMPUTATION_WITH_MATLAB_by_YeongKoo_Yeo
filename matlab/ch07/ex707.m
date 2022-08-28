% sldif.m
% x(1)=Ca2, x(2)=Ca3, ..., x(7)=Ca8
Dab = 1e-9; delx = 5e-4; ca10 = 1e-3; ca90 = 2e-3; ca0 = 6e-3; K = 1.5;
for k = 1:9, x0(k) = ca10 + (ca90 - ca10)*(k-1)/8; end
c0 = [x0(2) x0(3) x0(4) x0(5) x0(6) x0(7) x0(8)]; tf = 20000; tspan = [0 tf];
[t x] = ode45(@slf,tspan,c0,[],Dab,delx,ca0,ca10,ca90,K);
nt = length(t); xc1 = [ca10 ones(1,nt-1)*ca0/K]; % Ca1
xc9 = [ca90 ((4*x(2:end,7) - x(2:end,6))/3)']; % Ca9
xc = [];
for k = 1:7, xc = [xc x(:,k)]; end
c = [xc1' xc xc9'];
plot(t,c(:,3),t,c(:,5),':',t,c(:,7),'.-',t,c(:,9),'--'), xlabel('t(sec)'),
ylabel('C(kgmol/m^3)')
legend('C_{A3}','C_{A5}','C_{A7}','C_{A9}','location','best')

function dx = slf(t,x,Dab,delx,ca0,ca10,ca90,K)
% x(1)=Ca2, x(2)=Ca3, ..., x(7)=Ca8
if t == 0, ca1 = ca10; ca9 = ca90;
else, ca1 = ca0/K; ca9 = (4*x(7) - x(6))/3; end
dx = [Dab*(x(2)-2*x(1)+ca1)/(delx^2); Dab*(x(3)-2*x(2)+x(1))/(delx^2);
Dab*(x(4)-2*x(3)+x(2))/(delx^2); Dab*(x(5)-2*x(4)+x(3))/(delx^2);
Dab*(x(6)-2*x(5)+x(4))/(delx^2); Dab*(x(7)-2*x(6)+x(5))/(delx^2);
Dab*(ca9-2*x(7)+x(6))/(delx^2)];
end
