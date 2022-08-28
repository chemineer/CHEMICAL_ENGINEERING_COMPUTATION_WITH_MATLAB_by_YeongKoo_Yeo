function [u q]  =  hypbpde(f1,f2,g0,g1,xspan,tspan,nx,nt,alpa) 
% 1-dimensional hyperbolic PDE (wave equation: u_tt = alpa*u_xx) 
% Outputs: 
% u: solution matrix [u(x,t)] 
% Inputs: 
% f1, f2: initial conditions f1(x) and f2(x) (at t = 0) 
% g0, g1: boundary conditions g0(t) and g1(t) (at x = 0,1) 
% nx, nt: number of subintervals along x and t directions 
% xspan, tspan: range of x and t 
% alpa: coefficient 
% Output 
% u: solution matrix 
x0  =  xspan(1); xf  =  xspan(2); t0  =  tspan(1); tf  =  tspan(2); 
dx  =  (xf - x0)/nx; dt  =  (tf - t0)/nt; x  =  [0:nx]'*dx; t  =  [0:nt]*dt; 
q  =  alpa*(dt/dx)^2; q1  =  q/2; q2  =  2*(1-q); 
u(:,1)  =  f1(x); 
for k  =  1:nt+1, u([1 nx+1],k)  =  [g0(t(k)); g1(t(k))]; end 
u(2:nx,2) = q1*u(1:nx-1,1) + (1-q)*u(2:nx,1) + q1*u(3:nx+1,1) + dt*f2(x(2:nx)); 
for k  =  3:nt+1 
u(2:nx,k)  =  q*u(1:nx-1,k-1) + q2*u(2:nx,k-1) + q*u(3:nx+1,k-1) - u(2:nx,k-2); 
end 
surf(t,x,u) % display results 
colormap(gray); xlabel('t'), ylabel('x'), zlabel('u(t,x)') 
end 