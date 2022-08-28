% Solution of an ODE
tspan  =  [0 1]; y0  =  -1; dy  =  @(t,y) exp(-t); % define ODE 
[t y]  =  ode45(dy, tspan, y0); plot(t,y), grid, xlabel('t'), ylabel('y(t)') 