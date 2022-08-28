function x = congrad(A,x0,b) 
% Solves Ax=b by conjugate gradient method (A is square) 
% input: 
% A: coefficient matrix 
% b: constant vector (A*x=b) 
% x0: initial guess 
% Initialization 
itmax = 500; tol = 1e-6;  
b = (b(:).')'; x0 = (x0(:).')'; % b and x0 should be column vector 
if det(A) == 0   
fprintf('\n Rank = %7.3g\n',rank(A)); error('Matrix A is singular.'); 
end 
n = length(b); x = x0; M = A;  
for k = 1:n % check diagonal dominancy of A    
if sum(M(k,:)) > A(k,k), disp('A is not diagonally dominant.'); return; end 
end;  
% Iteration 
s = b - A*x; d = s; 
for k = 1:n    
v = A*d; alpa = dot(d,s)/dot(d,v); x = x + alpa*d; s = b - A*x;    
if sqrt(dot(s,s)) < tol, return;     
else, beta = -dot(s,v)/dot(d,v); d = s + beta*d; end 
end 
end 