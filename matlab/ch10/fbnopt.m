function [x, f, fint] = fbnopt(objfun, a, b, n) 
% fbnopt.m: 1-dimensional Fibonacci search method 
% inputs: 
%   objfun: objective function 
%   a, b: initial interval points 
% n: number of reduction 
% output: 
%   fint: size of final interval 
%   x: optimum point 
%   f: function value at the optimum point 
% Sample run: 
% a = 0; b = 50; n = 20; fobj= @(x) -870*x + 102*x^2 - 5*x^3; 
% [x, f, fint] = fbnopt(fobj, a, b, n); 

v = (sqrt(5)-1)/2; w = (1-sqrt(5))/(1+sqrt(5)); x1 = a; x4 = b; alpha = v*(1- w^n)/(1-w^(n+1)); 
x3 = alpha*x4 + (1-alpha)*x1; f3 = objfun(x3); 
for k = 1: n - 1    
if k == n - 1, x2 = 0.01*x1 + 0.99*x3;    
else, x2 = alpha*x1 + (1-alpha)*x4; end    
f2 = objfun(x2);    
if (f2 < f3), x4 = x3; x3 = x2; f3 = f2;     
else, x1 = x4; x4 = x2; f4 = f2; end    
alpha = v*(1- w^(n - k)) / (1-w^(n - k + 1)); 
end 
x = x3; f = f3; fint = abs(x1 - x4); 
end 