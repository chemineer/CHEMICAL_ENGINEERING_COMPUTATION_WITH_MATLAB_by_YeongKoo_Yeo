function [x,t,u] = parab1D(nx,nt,dx,dt,alpha,u0,bc,func,varargin)
% Solve 1-dimensional parabolic PDE
% Outputs:
% x, t; vectors of x and t values
% u: matrix of dependent variables [u(x,t)]
% Inputs:
% nx, nt: number of divisions in x and t direction
% dx, dy: x and t increment
% alpha: coefficient of equation
% u0: vector of u distribution at t=0
% bc: a matrix containing types and values of boundary conditions
% in x direction (2x2 or 2x3 matrix)
% row 1: conditions at lower x, row 2: conditions at upper x
% 1st column: type of condition
% (1: Dirichlet condition, values of u in 2nd column
% 2: Neumann condition, values of u' in 2nd column
% 3: Robbins condition, 2-3 columns contain constant and coef. of u)
% Initialization
if nargin < 7, error('Insufficient number of inputs.'); end
nx = fix(nx); x = [0:nx]*dx; nt = fix(nt); t = [0:nt]*dt; r = alpha*dt/dx^2;
u0 = (u0(:).')'; % confirm column vector
if length(u0) ~= nx+1, error('Incorrect length of the initial condition vector.'); end
[a,b] = size(bc);
if a ~= 2, error('Invalid number of boundary conditions.'); end
if b < 2 | b > 3, error('Invalid boundary condition.'); end
if b == 2 & max(bc(:,1)) <= 2, bc = [bc zeros(2,1)]; end
u(:,1) = u0; c = zeros(nx+1,1);
% Iteration according to t
for n = 2:nt+1 % lower x boundary condition
switch bc(1,1)
case 1, A(1,1) = 1; c(1) = bc(1,2);
case {2, 3}
A(1,1) = -3/(2*dx) - bc(1,3); A(1,2) = 2/dx; A(1,3) = -1/(2*dx); c(1) = bc(1,2);
end
% Interior points
for i = 2:nx
A(i,i-1) = -r; A(i,i) = 2*(1+r); A(i,i+1) = -r;
c(i) = r*u(i-1,n-1) + 2*(1-r)*u(i,n-1) + r*u(i+1,n-1);
if nargin >= 8 % Nonhomogeneous equation
intercept = feval(func,0,x(i),t(n),varargin{:});
slope = feval(func,1,x(i),t(n),varargin{:}) - intercept; A(i,i) = A(i,i) - dt*slope;
c(i) = c(i) + dt*feval(func,u(i,n-1),x(i),t(n-1),varargin{:}) + dt*intercept;
end
end
switch bc(2,1) % upper x boundary condition
case 1, A(nx+1,nx+1) = 1; c(nx+1) = bc(2,2);
case {2, 3}
A(nx+1,nx+1) = 3/(2*dx) - bc(2,3);
A(nx+1,nx) = -2/dx; A(nx+1,nx-1) = 1/(2*dx); c(nx+1) = bc(2,2);
end
u(:,n) = inv(A)*c;
end
end
