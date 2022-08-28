% pinelemtemp.m: temperature profiles of pin elements
Ta = 25; L0 = 0.01; D = 0.002; L = L0/5; Ac = pi*D*L; Ak = 3.14e-6; hc = 20; % data
rho = 2707; k = 220; Cp = 896; % properties of pure aluminum
V = Ak*L; Ch = rho*Cp*V; Rc = 1/(hc*Ac); Rce = 1/(hc*Ak); Rk = L/(k*Ak); Rk0 = L/(2*k*Ak);
T0 = Ta*ones(1,5); Tb = 100; % initial pin and base temperature at t = 0
tspan = [0 3];
dT = @(t,T) [(-(1/Rk0+1/Rk+1/Rc)*T(1)+T(2)/Rk+Ta/Rc+Tb/Rk0)/Ch;
(T(1)/Rk-(2/Rk+1/Rc)*T(2)+T(3)/Rk+Ta/Rc)/Ch;
(T(2)/Rk-(2/Rk+1/Rc)*T(3)+T(4)/Rk+Ta/Rc)/Ch;
(T(3)/Rk-(2/Rk+1/Rc)*T(4)+T(5)/Rk+Ta/Rc)/Ch;
(T(4)/Rk-(1/Rk+1/Rc+1/Rce)*T(5)+(1/Rc+1/Rce)*Ta)/Ch];
[t,T] = ode45(dT,tspan,T0);
T1 = T(:,1); T2 = T(:,2); T3 = T(:,3); T4 = T(:,4); T5 = T(:,5);
plot(t,T1,t,T2,':',t,T3,'.-',t,T4,'--',t,T5,'.'), xlabel('t(s)'), ylabel('T(deg.C)'), grid
legend('T_1','T_2','T_3','T_4','T_5')
% Steady-state temperature
Tst = @(x) [-(1/Rk0+1/Rk+1/Rc)*x(1)+x(2)/Rk+Ta/Rc+Tb/Rk0;
x(1)/Rk-(2/Rk+1/Rc)*x(2)+x(3)/Rk+Ta/Rc;
x(2)/Rk-(2/Rk+1/Rc)*x(3)+x(4)/Rk+Ta/Rc;
x(3)/Rk-(2/Rk+1/Rc)*x(4)+x(5)/Rk+Ta/Rc;
x(4)/Rk-(1/Rk+1/Rc+1/Rce)*x(5)+(1/Rc+1/Rce)*Ta];
x0 = Ta*ones(1,5); Ts = fsolve(Tst,x0);
fprintf('Steady-state temperature:\nPin element 1 = %g deg.C\n',Ts(1));
fprintf('Pin element 2 = %g deg.C\n',Ts(2)); fprintf('Pin element 3 = %g deg.C\n',Ts(3));
fprintf('Pin element 4 = %g deg.C\n',Ts(4)); fprintf('Pin element 5 = %g deg.C\n',Ts(5));
