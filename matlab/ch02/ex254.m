% flupb.m: fluidized packed bed catalytic reactor 
% IDASolve failed
A  =  0.17142; C  =  205.74; F  =  8000; Hg  =  320; Ht  =  266.67; Hw  =  1.6; % data 
Pe  =  0.1; Te  =  600; Tw  =  720; % data and parameters 
% define the differential equation system as an anonymous function 
dy = @(t,y) [Pe - y(1) + Hg*(y(3) - y(1)); Te - y(2) + Ht*(y(4) - y(2)) + Hw*(Tw - y(2)); Hg*(y(1) - y(3)*(1 + (6e-4 * exp(20.7 -15000/y(4)))))/A; Ht*((y(2) - y(4)) + F*(6e-4 * exp(20.7 -15000/y(4)))*y(3))/C]; 
tspan  =  [0 1500]; y0  =  [0.1 600 0 761]; % time range and initial conditions 
[t y]  =  ode15s(dy,tspan,y0); % use the stiff system solver ode15s 
subplot(1,2,1),plot(t,y(:,1),t,y(:,3),':'),grid,xlabel('\tau'),ylabel('Pressure(atm)'),legend('P','P_p') 
subplot(1,2,2),plot(t,y(:,2),t,y(:,4),':'),grid,xlabel('\tau'),ylabel('Temperature(R)'),legend('T','T_p') 