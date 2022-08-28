% ser3v.m 
v0 = 6; v = 12; V = 200; Ca0 = 2; Cb0 = 2; k = 0.5; tspan = 0:0.01:20; C0 = [Ca0 Cb0 0 0 0 0]; 
dCdt = @(t,C) [(v0*Ca0-v*C(1)-k*V*C(1)*C(2))/V; (v0*Cb0-v*C(2)-k*V*C(1)* C(2))/V;              
(v*C(1)-v*C(3)-k*V*C(3)*C(4))/V; (v*C(2)-v*C(4)-k*V*C(3)*C(4))/V;              
(v*C(3)-v*C(5)-k*V*C(5)*C(6))/V; (v*C(4)-v*C(6)-k*V*C(5)*C(6))/V]; 
[t,C]=ode45(dCdt, tspan, C0); plot(t,C(:,1),'-',t,C(:,3),':',t,C(:,5),'--'); 
xlabel('Time(min)'), ylabel('Concentration(gmol/dm^3)'), legend('C_ {A1}','C_{A2}','C_{A3}'); 