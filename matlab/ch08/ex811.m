% tempPDE.m
% Solve Laplace/Poisson equations using finite difference method
clear all;
distx = 1; % length of the plate in x direction (m)
disty = 1; % width of the plate in y direction (m)
ndx = 20; % number of divisions in x direction
ndy = 20; % number of divisions in y direction
rhf = -100e3/16; % right-hand side of the equation
% Boundary conditions: 1) Dirichlet, 2) Neumann, 3) Robbins
bc(1,1) = 3;  %  Lower x boundary condition: Robbins
bc(1,2) = -5*25; % constant (beta)
bc(1,3) = 5; % coefficient (gamma)
bc(2,1) = 3; % Upper x boundary condition: Robbins
bc(2,2) = 5*25; % constant (beta)
bc(2,3) = -5; % coefficient (gamma)
bc(3,1) = 3;  % Lower y boundary condition: Robbins
bc(3,2) = -5*25; % constant (beta)
bc(3,3) = 5; % coefficient (gamma)
bc(4,1) = 3; % Upper y boundary condition
bc(4,2) = 5*25; % constant (beta
bc(4,3) = -5; % coefficient (gamma)
[x,y,T] = ellipticPDE(ndx,ndy,distx/ndx,disty/ndy,bc,rhf); colormap(jet)
surf(y,x,T), xlabel('x(m)'), ylabel('y(m)'), zlabel('T(deg C)'), colorbar,
view(135,45)
