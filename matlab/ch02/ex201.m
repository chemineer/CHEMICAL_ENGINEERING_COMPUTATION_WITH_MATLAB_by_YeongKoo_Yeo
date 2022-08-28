% trilin.m: solve a linear system  
% Gauss elimination, Gauss-Seidel and conjugate gradient methods are used. 
A = zeros(6,6); for i = 1:6, A(i,i) = 2; end, for i = 1:5, A(i,i+1) = -1; A(i+1,i) = -1; end 
b = zeros(6,1); b(6,1) = 4; rho = 1; x0 = zeros(1,6); 
xg = gausselm(A,b); xs = GaussSeidel(A,b,rho); xc = congrad(A,x0,b); 
display('x by Gauss elimination method: '), display(xg'); 
display('x by Gauss-Seidel method: '), display(xs'); 
display('x by conjugate gradient method: '), display(xc'); 