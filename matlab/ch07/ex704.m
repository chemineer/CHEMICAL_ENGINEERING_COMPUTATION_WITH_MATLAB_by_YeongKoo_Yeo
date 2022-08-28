% drugf.m
x10 = 0.5; x20 = 0.4; x30 = 0.35; x40 = 0; x50 = 0;
V = 1200; tau = 240; Sa = 1; Sd = 0.4; rho = 1414.7;
tspan = [0 150]; x0 = [x10 x20 x30 x40 x50];
[t x] = ode45(@ruf,tspan,x0,[],V,rho,tau,Sa,Sd);
D1 = x(:,1); D2 = x(:,2); D3 = x(:,3); Cas = x(:,4); Cds = x(:,5);
subplot(1,2,1), plot(t,D1,t,D2,':',t,D3,'.-'), xlabel('t(min)'), ylabel('D(cm)')
legend('D_1','D_2','D_3'), axis tight
subplot(1,2,2), plot(t,Cas,t,Cds,':'), xlabel('t(min)'), ylabel('C(mg/cm^3)')
axis tight, legend('C_{AS}','C_{DS}','Location','best')
function dx = ruf(t,x,V,rho,tau,Sa,Sd)
% x(i)=D(i) (i=1,2,3), X(4)=C_AS, x(5)=C_Ds
sum0 = 0; sum1 = 0;
for i = 1:3
kL(i) = 1.2/x(i); Sw(i) = 0;
if x(i) > 0.3, fd(i) = -2*kL(i)*(Sa - x(4))/rho; Sw(i) = 1;
elseif x(i) <= 0.3 && x(i) >= 1e-5, fd(i) = -2*kL(i)*(Sd - x(5))/rho;
else, fd(i) = 0; end
sum0 = sum0 + Sw(i)*kL(i)*x(i); sum1 = sum1 + (1-Sw(i))*kL(i)*x(i);
end
fd(4) = pi*(Sa - x(4))*sum0/V - x(4)/tau; fd(5) = pi*(Sd - x(5))*sum1/V - x(4)/tau;
dx = [fd(1) fd(2) fd(3) fd(4) fd(5)]';
end
