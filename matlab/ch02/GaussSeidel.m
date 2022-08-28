function x = GaussSeidel(A,b,rho) 
% Solves Ax=b by Gauss-Seidel iterative method (A should be square and diagonally dominant) 
% input: 
% A: coefficient matrix 
% b: constant vector 
% rho: relaxation factor  

% Initialization 
itmax = 500; tol = 1e-8;  
b = (b(:).')'; % b should be column vector 
n = length(b); [nr nc] = size(A); 
if nr ~= nc, error('Matrix A is not square.'); end 
if nr ~= n, error('Matrix A and vector b are not consistent.'); end 
if det(A) == 0, fprintf('\n Rank = %7.3g\n',rank(A)); error('A is singular.'); end 
% Iteration 
M = A;  
for k = 1:n, M(k,k) = 0; x(k) = 0; end % initial values 
for k = 1:n % check diagonal dominancy of A    
if sum(M(k,:)) > A(k,k), disp('A is not diagonally dominant.'); return; end 
end; 
for k = 1:n, M(k,:) = M(k,:)/A(k,k); s(k) = b(k)/A(k,k); end 
iter = 0; x = (x(:).')'; 
while (1)    
x0 = x;     
for k = 1:n        
x(k) = rho*(s(k) - M(k,:)*x) + (1-rho)*x(k);        
if x(k) ~= 0, err(k) = abs((x(k) - x0(k))/x(k)); end    
end    
iter = iter + 1;    
if max(err) <= tol | iter >= itmax, break; end 
end 
end 