% rxncon.m 
clear all; 
k = 2.2; v0 = 0.05; V0 = 5; Cb0 = 0.025; Ca0 = 0.05; % data 
x0 = [Ca0 0 0 0]; tspan = [0 500]; % initial guess and time interval 
[t x] = ode45(@semibrx,tspan,x0,[],k,v0,V0,Cb0); % solve ODE system 
Ca = x(:,1); Cb = x(:,2); Cc = x(:,3); Cd = x(:,4); 
V = V0 + v0*t; Xa = (Ca0*V0 - Ca.*V)/(Ca0*V0); rA = k*Ca.*Cb; 
subplot(1,2,1), plot(t,Ca, t,Cb,':', t,Cc,'.-', t,Cd,'--') 
xlabel('t(sec)'), ylabel('Concentration(mol/dm^3)'), legend('C_A','C_B', 'C_C','C_D') 
subplot(1,2,2), plot(t,rA), xlabel('t(sec)'), ylabel('Reaction rate(mol/ dm^3sec)')  

function dxdt = semibrx(t,x,k,v0,V0,Cb0) 
% x(1)=Ca, x(2)=Cb, x(3)=Cc, x(4)=Cd 
rA = -k*x(1)*x(2); V = V0 + v0*t; 
dxdt = [rA - v0*x(1)/V; rA + (Cb0 - x(2))*v0/V; -rA - v0*x(3)/V; -rA - v0*x(4)/V]; 
end 