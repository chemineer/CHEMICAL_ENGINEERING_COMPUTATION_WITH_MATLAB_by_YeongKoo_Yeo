% van der Pol Equation
mu = 1; tspan = [0 25]; y0 = [1 1]; dy = @(t,y) [y(2); -y(1) + mu*(1-y(1)^2)*y(2)]; % define ODE 
[t y] = ode45(dy, tspan, y0); plot(t,y(:,1),t,y(:,2),':'), legend('y_1','y_2'), xlabel('t'), ylabel('y')