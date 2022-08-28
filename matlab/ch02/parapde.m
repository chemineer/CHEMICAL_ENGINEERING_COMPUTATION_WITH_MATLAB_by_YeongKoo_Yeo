function u  =  parapde(f,g0,g1,tf,nx,nt,alpa) 
% Solution of parabolib PDE (u_t = u_xx) using explicit method 
% input 
% f: initial condition u(x,0) = f(x) 
% g0,g1: boundary conditions u(0,t) = g0(t), u(1,t) = g1(t) 
% tf: final time 
% nx,nt: number of subintervals 
% Output 
% u: solution matrix 
h  =  1/nx; d  =  tf/nt; r  =  alpa*d/h^2; % r should be less than 1/2 
x  =  0:h:1; t  =  0:d:tf; 
u(1:nx+1,1)  =  f(x)'; % initial condition: u(x,0) = f(x) 
u(1,1:nt+1)  =  g0(t); u(nx+1,1:nt+1)  =  g1(t); % boundary conditions 
for k  =  1:nt 
u(2:nx,k+1)  =  r*u(1:nx-1,k) + (1-2*r)*u(2:nx,k) + r*u(3:nx+1,k); 
end 
u  =  u'; surf(x,t,u) % display results (x from left to right) 
colormap(gray); xlabel('x'), ylabel('t'), zlabel('T(t,x)') 
end 