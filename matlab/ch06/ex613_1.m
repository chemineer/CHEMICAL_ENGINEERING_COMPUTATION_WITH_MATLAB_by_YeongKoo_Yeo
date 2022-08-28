% vandrxn.m  
V = 1; F = 25; Caf = 10; k1 = 50; k2 = 100; k3 = 10; % data 
df = @(t,C) [-k1*C(1) - k3*C(1)^2 + F*(Caf - C(1))/V; k1*C(1) - k2*C(2) - F*C(2)/V]; 
tspan = [0 0.06]; C0 = [10 0]; [t C] = ode45(df,tspan,C0); % solve ODEs by ode45 
plot(t,C(:,1),t,C(:,2),':'), legend('C_A','C_B'), xlabel('t'), ylabel('C(t)') 