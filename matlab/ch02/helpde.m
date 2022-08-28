function [u,x,y]  =  helpde(f,g,qx0,qxf,qy0,qyf,R,m,n,crit,kmax) 
% Solve the Helmholtz equation: u_xx + u_yy + g(x,y)u(x,y)  =  f(x,y) 
% over the region R  =  [x0 xf y0 yf] 
% Boundary conditions: 
% u(x0,y)  =  qx0(y), u(xf,y)  =  qxf(y), u(x,y0)  =  qy0(x), u(x,yf)  =  qyf(x) 
% m: number of subintervals along x-axis 
% n: number of subintervals along y-axis 
% crit: tolerance 
% kmax: the maximum number of iterations 
x0  =  R(1); xf  =  R(2); y0  =  R(3); yf  =  R(4); 
dx  =  (xf - x0)/m; x  =  x0 + [0:m]*dx; 
dy  =  (yf - y0)/n; y  =  y0 + [0:n]*dy; 
dx2  =  dx*dx; dy2  =  dy*dy; dxy2  =  2*(dx2 + dy2); 
bx  =  dx2/dxy2; by  =  dy2/dxy2; bxy  =  bx*dy2; 
m1  =  m + 1; n1  =  n + 1; 
% Boundary conditions 
for k  =  1:n1, u(k,[1 m1])  =  [qx0(y(k)) qxf(y(k))]; end % at x = x0 and x = xf 
for k  =  1:m1, u([1 n1], k)  =  [qy0(x(k)) qyf(x(k))]; end % at y = y0 and y = yf 
% Initialization 
u(2:m,2:n)  =  0; 
for i  =  1:m 
for j  =  1:n, fv(i,j)  =  f(x(i),y(j)); gv(i,j)  =  g(x(i),y(j)); end 
end 
% Solve the difference equation by successive substitution method 
for k  =  1:kmax % perform iterative calculation 
for i  =  2:m 
for j  =  2:n 
u(i,j) =  by*(u(i+1,j)+u(i-1,j)) + bx*(u(i,j+1)+u(i,j-1)) + bxy*(gv(i,j)*u(i,j)-fv(i,j)); 
end 
end 
if k>1 & max(max(abs(u-u0)))<crit, break; end 
u0  =  u; % back substitution 
end 
end 