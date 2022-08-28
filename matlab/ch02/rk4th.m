function [t,y]  =  rk4th(f,tspan,y0,n) 
% Implements 4th-order Runge-Kutta method to solve dy/dt  =  f(t,y) 
% input 
% f: dy/dt 
% tspan: vector of initial and final values of independent variable 
% y0: initial value of dependent variable 
% n: number of subintervals 
% Output 
% t: vector of independent variable 
% y: vector of dependent variable (solution vector) 
t0  =  tspan(1); tf  =  tspan(2); h  =  (tf-t0)/n; 
t  =  (t0:h:tf)'; nt  =  length(t); y  =  y0*ones(nt,1); 
for k  =  1:nt-1 
k1  =  f(t(k),y(k)); k2  =  f(t(k)+h/2, y(k)+h*k1/2); 
k3  =  f(t(k)+h/2,y(k)+h*k2/2); k4  =  f(t(k)+h, y(k)+h*k3); 
y(k+1)  =  y(k) + h*(k1 + 2*k2 + 2*k3 + k4)/6; 
end 
end 