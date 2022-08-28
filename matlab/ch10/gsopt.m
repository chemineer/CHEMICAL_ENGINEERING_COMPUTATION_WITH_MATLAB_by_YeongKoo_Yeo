function [x, f, n] = gsopt(objfun, x1, h, crit)
% gsopt.m: golden section search method
% input:
% objfun: objective function
% x1: initial point
% h: step size
% crit: stopping criterion
% output:
% x: optimal point
% f: function value at the optimal point
% n: number of function evaluations
% sample run:
% crit=1e-6;,x1=1; h=0.1;
% fun = @(x) 85*(1-x^3)^2 + (1-x^2) + 3*(1-x)^2;
% [x, f, n] = gsopt(fun, x1, h, crit)

% Initialization
tau = (sqrt(5) - 1)/2; n = 0; f1 = objfun(x1); n = n+1;
x2 = x1 + h; f2 = objfun(x2); n = n+1;
% Golden section search
if f2 > f1, temp = x1; x1 = x2; x2 = temp; temp = f1; f1 = f2; f2 = temp; h = -h; end
while x1 < 1e50
h = h/tau; x4 = x2 + h; f4 = objfun(x4); n = n+1;
if f4 > f2, break; end
f1 = f2; x1 = x2; f2 = f4; x2 = x4;
end
fold = (f1 + f2 + f4) / 3; ind = 0;
while x1 < 1e50
if abs(x4-x1) < crit, break; end
x3 = tau*x4 + (1-tau)*x1; f3 = objfun(x3); n = n+1;
if f2 < f3, x4 = x1; x1 = x3; f4 = f1; f1 = f3;
else, x1 = x2; x2 = x3; f1 = f2; f2 = f3; end
fpr = (f1 + f2 + f4) / 3;
if abs(fpr - fold) < crit, ind = ind + 1; if ind == 2, break; end
else, ind = 0; end
fold = fpr;
end
x = x2; f = f2;
end
