function [u r] = parabDbc(nx,nt,dx,dt,alpha,u0,bci,bcf)
% Solve 1-dimensional parabolic PDE (constant Dirichlet boundary conditions)
% Outputs:
% u: matrix of dependent variables [u(x,t)]
% Inputs:
% nx, nt: number of divisions in x and t direction
% dx, dy: x and t increments
% alpha: coefficient of equation
% u0: row vector of u distribution at t=0
% bci, bcf: boundary conditions at both ends of x (scalar)

r = alpha*dt/dx^2; A = zeros(nx-1,nx-1); u = zeros(nt+1,nx+1); u(:,1) = bci*ones(nt+1,1);
u(:,nx+1) = bcf*ones(nt+1,1); u(1,:) = u0; A(1,1) = 1 + 2*r; A(1,2) = -r;
for i = 2:nx-2, A(i,i) = 1 + 2*r; A(i,i-1) = -r; A(i,i+1) = -r; end
A(nx-1,nx-2) = -r; A(nx-1,nx-1) = 1 + 2*r; b(1,1) = u0(2) + u0(1)*r;
for i = 2:nx-2, b(i,1) = u0(i+1); end
b(nx-1,1) = u0(nx) + u0(nx+1)*r; [L U] = lu(A);
for j = 2:nt+1
y = L\b; x = U\y; u(j,2:nx) = x'; b = x;
b(1,1) = b(1,1) + bci*r; b(nx-1,1) = b(nx-1,1) + bcf*r;
end
end

