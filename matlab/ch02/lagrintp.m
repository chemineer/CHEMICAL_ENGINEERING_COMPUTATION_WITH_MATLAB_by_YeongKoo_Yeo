function yi  =  lagrintp(x,y,xi) 
% Implements Lagrange interpolation method 
% input 
% x: vector of independent variables 
% y: vector of dependent variables 
% xi: vector of independent variable at which the interpolation is performed 
% Output 
% yi: vector of interpolated values at xi 
% Determine coefficients of Lagrange polynomial 
n  =  length(x); m  =  length(xi); 
if n ~=  length(y), error('x and y should be the same length.'); end 
for k  =  1:n 
dx(k)  =  1; 
for j  =  1:n 
if j ~=  k, dx(k)  =  dx(k)*(x(k) - x(j)); end 
a(k)  =  y(k)/dx(k); 
end 
end 
% Interpolation using Lagrange polynomial 
for i  =  1:m, yi(i)  =  0; 
for j  =  1:n, q(j)  =  1; 
for k  =  1:n, if (j ~=  k), q(j)  =  q(j)*(xi(i) - x(k)); end; end 
yi(i)  =  yi(i) + a(j)*q(j); 
end 
end 
end 