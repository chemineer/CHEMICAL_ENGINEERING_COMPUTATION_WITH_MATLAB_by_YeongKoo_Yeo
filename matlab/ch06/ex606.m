% serrb.m 
v = 25; k1 = 0.2; k2 = 0.1; Ca0 = 2.5; Cb0 = 0; Cc0 = 0; tspan = 0:0.01:15; C0 = [Ca0 Cb0 Cc0]; 
dCdt = @(t,C) [-k1*C(1); k1*C(1)-k2*C(2); k2*C(2)]; [t,C]=ode45(dCdt, tspan, C0); 
plot(t,C(:,1),'-',t,C(:,2),':',t,C(:,3),'--');  
xlabel('Time(min)'), ylabel('Concentration(mol/liter)'), legend('A','B','C'); 
Cmax = max(C(:,2)), tmax = t(find(C(:,2) == Cmax)), Vol = v*tmax 