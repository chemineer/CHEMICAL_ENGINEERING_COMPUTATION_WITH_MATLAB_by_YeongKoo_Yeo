% fembat.m: batch bioreactor for ethanol production
F = 1; Sf = 100; % data
dx = @(t,x) [F;
(0.408*x(3)*exp(-0.028*x(4))/(0.22+x(3)) - F/x(1))*x(2);
-4.08*x(3)*x(2)*exp(-0.028*x(4))/(0.22+x(3)) + F*(Sf - x(3))/x(1);
x(3)*x(2)*exp(-0.015*x(4))/(0.44+x(3)) - F*x(4)/x(1)];
tint = [0 20]; x0 = [1 0.2 100 0];
[t X] = ode45(dx,tint,x0); % x(1)=V, x(2)=x, x(3)=S, x(4)=P
x = X(:,2); S = X(:,3); P = X(:,4);
mu = 0.408*S.*exp(-0.028*P)./(0.22 + S);
phi = S.*exp(-0.015*P)./(0.44 + S);
subplot(1,2,1), plot(t,x,t,S,'--',t,P,':'), xlabel('t(hr)')
ylabel('Concentration(g/liter)'), grid, legend('x','S','P','location','best')
subplot(1,2,2), plot(t,mu,t,phi,'--'), xlabel('t(hr)')
ylabel('\mu and \pi'), grid, legend('\mu','\pi','location','best')
