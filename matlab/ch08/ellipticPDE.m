function [x,y,U] = ellipticPDE(nx,ny,dx,dy,bc,f)
% Solve 2-dimensional elliptic PDE
% output:
% x, y: vectors of x and y values
% U: matrix of dependent variables (U(x,y))
% input:
% nx, ny: number of divisions in x and y direction
% dx, dy: x and y increments
% f: constant for Poisson equation
% bc: a matrix containing types and values of boundary conditions
% in x and y directions (4x2 or 4x3 matrix)
% order of appearing: lower x, upper x, lower y, upper y
% in rows 1 to 4 of the matric bc
% 1st column: type of condition:
% (1: Dirichlet condition, values of u in 2nd column
% 2: Neumann condition, values of u' in 2nd column
% 3: Robbins condition, 2-3 columns contain constant and coef. of u)

% Initialization
if nargin < 5, error('Invalid number of inputs.'); end
[a,b] = size(bc);
if a ~= 4, error('Invalid number of boundary conditions.'); end
if b < 2 | b > 3, error('Invalid boundary condition.'); end
if b == 2 & max(bc(:,1)) <= 2, bc = [bc zeros(4,1)]; end
if nargin < 6 | isempty(f), f = 0; end
nx = fix(nx); x = [0:nx]*dx; ny = fix(ny); y = [0:ny]*dy; dx2 = 1/dx^2; dy2 = 1/dy^2;
% Coefficient matrix and constant vector
n = (nx+1)*(ny+1); A = zeros(n); c = zeros(n,1);
onex = diag(diag(ones(nx-1))); oney = diag(diag(ones(ny-1)));
% Interior nodes
i = [2:nx];
for j = 2:ny
ind = (j-1)*(nx+1)+i; A(ind,ind) = -2*(dx2+dy2)*onex;
A(ind,ind+1) = A(ind,ind+1) + dx2*onex; A(ind,ind-1) = A(ind,ind-1) + dx2*onex;
A(ind,ind+nx+1) = A(ind,ind+nx+1) + dy2*onex;
A(ind,ind-nx-1) = A(ind,ind-nx-1) + dy2*onex; c(ind) = f*ones(nx-1,1);
end
% Lower x boundary condition
switch bc(1,1)
case 1
ind = ([2:ny]-1)*(nx+1)+1; A(ind,ind) = A(ind,ind) + oney; c(ind) = bc(1,2)*ones(ny-1,1);
case {2, 3}
ind = ([2:ny]-1)*(nx+1)+1; A(ind,ind) = A(ind,ind) - (3/(2*dx) + bc(1,3))*oney;
A(ind,ind+1) = A(ind,ind+1) + 2/dx*oney; A(ind,ind+2) = A(ind,ind+2) - 1/ (2*dx)*oney;
c(ind) = bc(1,2)*ones(ny-1,1);
end
% Upper x boundary condition
switch bc(2,1)
case 1
ind = [2:ny]*(nx+1); A(ind,ind) = A(ind,ind) + oney; c(ind) = bc(2,2)*ones(ny-1,1);
case {2, 3}
ind = [2:ny]*(nx+1); A(ind,ind) = A(ind,ind) + (3/(2*dx) - bc(2,3))*oney;
A(ind,ind-1) = A(ind,ind-1) - 2/dx*oney;
A(ind,ind-2) = A(ind,ind-2) + 1/(2*dx)*oney; c(ind) = bc(2,2)*ones(ny-1,1);
end
% Lower y boundary condition
switch bc(3,1)
case 1
ind = [2:nx]; A(ind,ind) = A(ind,ind) + onex; c(ind) = bc(3,2)*ones(nx-1,1);
case {2, 3}
ind = [2:nx]; A(ind,ind) = A(ind,ind) - (3/(2*dy) + bc(3,3))*onex;
A(ind,ind+nx+1) = 2/dy*onex; A(ind,ind+2*(nx+1)) = -1/(2*dy)*onex;
c(ind) = bc(3,2)*ones(nx-1,1);
end
% Upper y boundary condition
switch bc(4,1)
case 1
ind = ny*(nx+1)+[2:nx]; A(ind,ind) = A(ind,ind) + onex; c(ind) = bc(4,2)*ones(nx-1,1);
case {2, 3}
ind = ny*(nx+1)+[2:nx]; A(ind,ind) = A(ind,ind) + (3/(2*dy) - bc(4,3))*onex;
A(ind,ind-(nx+1)) = A(ind,ind-(nx+1)) - 2/dy*onex;
A(ind,ind-2*(nx+1)) = A(ind,ind-2*(nx+1)) + 1/(2*dy)*onex; c(ind) = bc(4,2)*ones(nx-1,1);
end
% Corner nodes
A(1,1) = 1; A(1,2) = -1/2; A(1,nx+2) = -1/2; c(1) = 0;
A(nx+1,nx+1) = 1; A(nx+1,nx) = -1/2; A(nx+1,2*(nx+1)) = -1/2; c(nx+1) = 0;
A(ny*(nx+1)+1,ny*(nx+1)+1) = 1; A(ny*(nx+1)+1,ny*(nx+1)+2) = -1/2;
A(ny*(nx+1)+1,(ny-1)*(nx+1)+1) = -1/2; c(ny*(nx+1)+1) = 0;
A(n,n) = 1; A(n,n-1) = -1/2; A(n,n-(nx+1)) = -1/2; c(n) = 0;
u = inv(A)*c; % solve equation
% Rearrange results in matrix form
for k = 1:ny+1, U(k,1:nx+1) = u((k-1)*(nx+1)+1:k*(nx+1))'; end
end
