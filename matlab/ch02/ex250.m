% Well-Mixed Tanks

tspan  =  [0 2]; y0  =  [0.15 0.15 0.15 20]; 
dy  =  @(t,y) [(2.25 -27*y(1))./y(4); 15*(y(1)-y(2))/20; 15*(y(2)-y(3))/20; 12]; % equation system 
[t y]  =  ode45(dy, tspan, y0); plot(t,y(:,1),t,y(:,2),'.-',t,y(:,3),':'), xlabel('t(min)'), ylabel('x(t)') 
legend('x_1', 'x_2', 'x_3', 'location', 'best') 