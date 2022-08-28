% Numerical Solution of ODE 
f  =  @(t,y) 5*exp(0.6*t) - 2*y; y0  =  1.5; tspan  =  [0 3]; n  =  5; 
[t,y]  =  eulerde(f,tspan,y0,n); disp([t y]) % explicit Euler method 
[t,y]  =  rk4th(f,tspan,y0,n); disp([t y]) % 4th-order Runge-Kutta method 