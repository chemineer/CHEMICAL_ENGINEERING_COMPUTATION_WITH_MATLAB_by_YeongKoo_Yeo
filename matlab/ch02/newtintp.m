function yi  =  newtintp(x,y,xi) 
% Implements Newton interpolation method 
% input 
% x: vector of independent variables 
% y: vector of dependent variables 
% xi: vector of independent variable at which the interpolation is performed 
% Output 
% yi: vector of interpolated values at xi 
% Determine coefficients of Newton polynomial 
n  =  length(x); m  =  length(xi); 
if n ~=  length(y), error('x and y should be the same length.'); end 
a(1)  =  y(1); 
for k  =  1:n-1, df(k,1)  =  (y(k+1) - y(k))/(x(k+1) - x(k)); end 
for j  =  2:n-1 
for k  =  1:n-j, df(k,j)  =  (df(k+1,j-1) - df(k,j-1))/(x(k+j) - x(k)); end 
end 
for k  =  2:n, a(k)  =  df(1,k-1); end 
% Interpolation using Newton polynomial 
for k  =  1:m 
s(1)  =  1; yi(k)  =  a(1); 
for j  =  2:n, s(j)  =  (xi(k) - x(j-1))*s(j-1); yi(k)  =  yi(k) + a(j)*s(j); end 
end 
end 