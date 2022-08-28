% penrxn.m: Penicillin production reaction
a  =  [13.2 0.95 1.76]; x0  =  [0.028 0]; tspan  =  [0 1]; 
dx  =  @(t,x) [a(1)*x(1) - (a(1)/a(2))*x(1)^2; a(3)*x(1)]; [t x]  =  ode45(dx,tspan,x0); 
% Generate profiles 
plot(t,x(:,1),t,x(:,2),'--'), grid, legend('x_1','x_2','location','best') 
xlabel('t'), ylabel('x_1 and x_2') 