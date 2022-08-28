function p  =  lspolfit(x,y,m) 
% Implements least-squares polynomial fitting 
% input 
% x: vector of independent variables 
% y: vector of measured dependent variables 
% m: order of polynomial 
% Output 
% p: coefficient vector of mth order polynomial 
n  =  length(x); m  =  m+1; 
if n ~=  length(y), error('x and y should be the same length.'); end 
A  =  zeros(m,m); 
% Define matrix A 
for k  =  1:2*(m-1), s(k)  =  sum(x.^k); end; 
A(1,1)  =  n; A(1,2:m)  =  s(1:m-1); 
for k  =  2:m, A(k,:)  =  s(k-1:k-2+m); end; 
% Define vector b 
for k  =  1:m, b(k)  =  sum((x.^(k-1)).*y); end; b  =  b'; 
% Determine parameter using backslash operator 
c  =  A\b; p  =  c(end:-1:1); % coefficients in descending order 
end 