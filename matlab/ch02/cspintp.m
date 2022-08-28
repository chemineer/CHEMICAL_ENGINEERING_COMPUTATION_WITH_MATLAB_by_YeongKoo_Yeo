function yi  =  cspintp(x,y,xi) 
% Implements cubic spline interpolation method 
% input 
% x: vector of independent variables 
% y: vector of dependent variables 
% xi: vector of independent variable at which the interpolation is performed 
% Output 
% yi: vector of interpolated values at xi 
% Formulation of tridiagonal system 
n  =  length(x); m  =  length(xi); 
if n ~=  length(y), error('x and y should be the same length.'); end 
for k  =  1:n-1 % parameters to form tridiagonal system 
h(k)  =  x(k+1) - x(k); w(k)  =  (y(k+1) - y(k))/h(k); 
end 
for k  =  1:n-2 % right-hand side and diagonal elements 
v(k)  =  w(k+1) - w(k); d(k)  =  2*(h(k) + h(k+1)); 
end 
for k  =  1:n-3, U(k)  =  h(k+1); L(k+1)  =  U(k); end 
L(1)  =  0; U(n-2)  =  0; 
% Solve tridiagonal system to find coefficients of spline functions 
v(1)  =  v(1)/d(1); U(1)  =  U(1)/d(1); 
for k  =  2:n-3 
dn  =  d(k)-L(k)*U(k-1); U(k)  =  U(k)/dn; v(k)  =  (v(k)-L(k)*v(k-1))/dn; 
end 
v(n-2) = (v(n-2) - L(n-2)*v(n-3))/(d(n-2) - L(n-2)*U(n-3)); a(n-2) = v(n-2); 
for k  =  n-3:-1:1, a(k)  =  v(k) - U(k)*a(k+1); end; % coefficient vector a 
b(1)  =  y(1)/h(1); c(1)  =  y(2)/h(1) - a(1)*h(1); % coefficient vectors b and c 
for k  =  2:n-2, b(k)  =  y(k)/h(k) - a(k-1)*h(k); c(k)  =  y(k+1)/h(k) - a(k)*h(k); end 
b(n-1)  =  y(n-1)/h(n-1) - a(n-2)*h(n-1); c(n-1)  =  y(n)/h(n-1); 
a  =  [0 a]; % a(1) should be 0 
% Interpolation at xi 
for k  =  1:m 
for j  =  1:n-1, if xi(k) >=  x(j), id  =  j; end; end 
if xi(k) > x(n), error('xi > max(x): cannot interpolate.'); end 
if xi(k) < x(1), error('xi < min(x): cannot interpolate.'); end 
h  =  x(id+1) - x(id); 
if id  ==  1 
s(k)  =  a(2)*(xi(k)-x(1))^3/h + b(1)*(x(2)-xi(k)) + c(1)*(xi(k)-x(1)); 
elseif id  ==  n-1 
s(k)  =  a(n-1)*(x(n)-xi(k))^3/h + b(n-1)*(x(n)-xi(k)) + c(n-1)*(xi(k)-x(n-1)); 
else 
s(k)  =  (a(id)*(x(id+1)-xi(k))^3 + a(id+1)*(xi(k)-x(id))^3)/h... 
+ b(id)*(x(id+1)-xi(k)) + c(id)*(xi(k)-x(id)); 
end 
end 
yi  =  s; % vector of interpolated values at xi 
end 