function z  =  trapzoid(f,a,b,n) 
% Implements trapezoidal rule to integrate f(x) from a to b 
% input 
% f: function handle to be integrated 
% a,b: limits of integration 
% n: number of subinterval points 
% Output 
% z: integral of f(x) 
h  =  (b-a)/n; s  =  f(a); 
for k  =  1:n-1, x  =  a + h*k; s  =  s + 2*f(x); end 
s  =  s + f(b); z  =  h*s/2; 
end 
function s  =  trapzoidat(x,y) 
% Implements trapezoidal rule to integrate data set (x,y) 
% input 
% x: vector of independent variables 
% y: vector of dependent variables 
% Output 
% s: integral of y 
n  =  length(x); s  =  0; 
for k  =  1:n-1, s  =  s + (y(k) + y(k+1))*(x(k+1) - x(k))/2; end 
end 