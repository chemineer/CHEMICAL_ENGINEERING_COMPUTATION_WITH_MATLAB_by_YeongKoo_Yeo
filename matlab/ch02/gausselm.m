function x = gausselm(A,b) 
% Solves Ax=b by Gauss elimination method (A: square(m=n)) 
b = (b(:).')';  % b must be a column vector 
n = length(b); 
% Elimination phase 
for k = 1:n-1    
for i = k+1:n        
if A(i,k) ~= 0            
c = A(i,k)/A(k,k); A(i,k+1:n) = A(i,k+1:n) - c*A(k,k+1:n);            
b(i) = b(i) - c*b(k);        
end    
end 
end 
% Solution phase: back substitution 
for k = n:-1:1, b(k) = (b(k) - A(k,k+1:n)*b(k+1:n))/A(k,k); end 
x = b; % vector b contains the solution 
end 