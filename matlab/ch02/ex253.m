% stiffode.m: solution of stiff differential equations 
tspan  =  [0 20]; y0  =  [0.05 5]; k  =  0.3; K  =  1e-6; % data 
dy  =  @(t,y) [k*y(1)*y(2)/(K+y(2)); -0.75*k*y(1)*y(2)/(K+y(2))]; % define equations 
[t y]  =  ode15s(dy, tspan, y0); % use the built-in solver ode15s 
plot(t,y(:,1),t,y(:,2),':'), grid, xlabel('t(min)'), ylabel('y(t)'), legend('B(t)','S(t)','location','best') 