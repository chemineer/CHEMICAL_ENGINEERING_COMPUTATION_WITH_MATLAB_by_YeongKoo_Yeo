function [x,y,t,u] = parab2D(nx,ny,nt,dx,dy,dt,alpha,u0,bc,func,varargin)
% Solve 2-dimensional parabolic PDE: use explicit formula
% Outputs:
% x, y, t: vectors of x, y and t values
% u: 3D array of dependent variables [u(x,y,t)]
% Inputs:
% nx, ny, nt: number of divisions in x, y and t direction
% dx, dy: x and y increments
% dt: t increments (leave empty to use the default values)
% alpha: coefficient of equation
% u0: matrix of u distribution at t=0 [u0(x,y)]
% bc: a matrix containing types and values of boundary conditions
% in x and y directions (4x2 or 4x3 matrix)
% order of appearing: lower x, upper x, lower y, upper y
% in rows 1 to 4 of the matric bc
% 1st column: type of condition
% (1: Dirichlet condition, values of u in 2nd column
% 2: Neumann condition, values of u' in 2nd column
% 3: Robbins condition, 2-3 columns contain constant and coef. of u)

% Initialization
if nargin < 9, error(' Invalid number of inputs.'); end
nx = fix(nx); x = [0:nx]*dx; ny = fix(ny); y = [0:ny]*dy;
% Check dt for stability
tmax = dt*nt;
if isempty(dt) | dt > (dx^2+dy^2)/(16*alpha)
dt = (dx^2+dy^2)/(16*alpha); nt = tmax/dt+1;
fprintf('\ndt is adjusted to %6.2e (nt=%3d)\n',dt,fix(nt))
end
nt = fix(nt); t = [0:nt]*dt; rx = alpha*dt/dx^2; ry = alpha*dt/dy^2;
[r0,c0] = size(u0);
if r0 ~= nx+1 | c0 ~= ny+1
error('Incorrect size of the initial condition matrix.')
end
[a,b] = size(bc);
if a ~= 4, error('Invalid number of boundary conditions.'); end
if b < 2 | b > 3, error('Invalid boundary condition.'); end
if b == 2 & max(bc(:,1)) <= 2, bc = [bc zeros(4,1)]; end
% Solution of PDE
u(:,:,1) = u0;
for n = 1:nt
for i = 2:nx
for j = 2:ny
u(i,j,n+1) = rx*(u(i+1,j,n)+u(i-1,j,n))+ ry*(u(i,j+1,n) + u(i,j-1,n)) + (1-2*rx-2*ry)*u(i,j,n);
if nargin >= 10
u(i,j,n+1) = u(i,j,n+1) + dt*feval(func,u(i,j,n),x(i),y(j),t(n),varargin{:});
end
end
end
% Lower x boundary condition
switch bc(1,1)
case 1, u(1,2:ny,n+1) = bc(1,2) * ones(1,ny-1,1);
case {2, 3}
u(1,2:ny,n+1) = (-2*bc(1,2)*dx + 4*u(2,2:ny,n+1) - u(3,2:ny,n+1)) / (2*bc(1,3)*dx + 3);
end
% Upper x boundary condition
switch bc(2,1)
case 1, u(nx+1,2:ny,n+1) = bc(2,2) * ones(1,ny-1,1);
case {2, 3}
u(nx+1,2:ny,n+1) = (-2*bc(2,2)*dx - 4*u(nx,2:ny,n+1) + u(nx-1,2:ny,n+1)) / (2*bc(2,3)*dx - 3);
end
% Lower y boundary condition
switch bc(3,1)
case 1, u(2:nx,1,n+1) = bc(3,2) * ones(nx-1,1,1);
case {2, 3}
u(2:nx,1,n+1) = (-2*bc(3,2)*dy + 4*u(2:nx,2,n+1) - u(2:nx,3,n+1)) / (2*bc(3,3)*dy + 3);
end
% Upper y boundary conditionCorner nodes
switch bc(4,1)
case 1, u(2:nx,ny+1,n+1) = bc(4,2) * ones(nx-1,1,1);
case {2, 3}
u(2:nx,ny+1,n+1) = (-2*bc(4,2)*dy - 4*u(2:nx,ny,n+1) + u(2:nx,ny-1,n+1)) / (2*bc(4,3)*dy - 3);
end
end
% Corner nodes
u(1,1,:) = (u(1,2,:) + u(2,1,:)) / 2;
u(nx+1,1,:) = (u(nx+1,2,:) + u(nx,1,:)) / 2;
u(1,ny+1,:) = (u(1,ny,:) + u(2,ny+1,:)) / 2;
u(nx+1,ny+1,:) = (u(nx+1,ny,:) + u(nx,ny+1,:)) / 2;
end
