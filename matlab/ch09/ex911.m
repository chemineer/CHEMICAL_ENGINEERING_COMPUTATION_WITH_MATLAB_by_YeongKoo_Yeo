% hybridtanks.m: hybrid two-tank heater system
% Data structure
ht.A1 = 0.25; ht.A2 = 0.25; ht.F0 = 0.4; ht.T0 = 25; ht.H = 0.5;
ht.c1 = 0.6; ht.c2 = 0.6; ht.rCp = 4180; ht.Q1 = 6000; ht.Q2 = 6000;
% Initial guesses
h10 = 0.4; h20 = 0.35; z0 = [h10; h20; 25; 25]; tspan = [0 10];
% Solve DE system
[t,z] = ode45(@LTmodel,tspan,z0,[],ht); h1 = z(:,1); h2 = z(:,2); T1 = z(:,3);
T2 = z(:,4);
% Plot results
subplot(1,2,1), plot(t,h1,t,h2,'--'),xlabel('t(min)'), ylabel('h_1,h_2(m)')
legend('h_1','h_2','location','best')
subplot(1,2,2), plot(t,T1,t,T2,'--'),xlabel('t(min)'), ylabel('T_1,T_2 (deg.C)'), legend('T_1','T_2','location','best')
