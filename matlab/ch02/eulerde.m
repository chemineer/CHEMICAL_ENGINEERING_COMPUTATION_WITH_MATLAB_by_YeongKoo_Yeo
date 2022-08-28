function [t,y]  =  eulerde(f,tspan,y0,n) 
% Implements explicit Euler method to solve dy/dt  =  f(t,y) 
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
for k  =  2:nt, y(k)  =  y(k-1) + h*f(t(k-1),y(k-1)); end 
end 