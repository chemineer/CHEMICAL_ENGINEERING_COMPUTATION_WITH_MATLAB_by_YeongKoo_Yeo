function z  =  newtrapmv(fun,x0) 
% Implements Newton-Raphson method to find roots of a system of equations 
% input 
% fun: function handle that returns f(x) = [f1,f2,...,fn] 
% x0: vector of initial guesses 
% Output 
% z: zeroes of f(x) 
tol  =  1e-8; kmax  =  1e3; 
x  =  x0; if size(x,1)  ==  1, x  =  x'; end % x should be column vector 
for k  =  1:kmax 
[J,f]  =  jacob(fun,x); 
if sqrt(dot(f,f)/length(x)) < tol, z  =  x; return; end 
dx  =  J\(-f); x  =  x + dx; 
if sqrt(dot(dx,dx)/length(x)) < tol*max(abs(x),1), z  =  x; return; end 
if k >= kmax, disp('The Newton-Raphson method does not converge.'); break; end 

end 
fprintf('Number of iterations: %g\n', k); 
end 
function [J,f]  =  jacob(fun,x) 
% Computes the Jacobian matrix J and f(x) 
hx  =  1e-4; n  =  length(x); J  =  zeros(n); f0  =  feval(fun,x); f  =  f0; 
for k  =  1:n 
tempx  =  x(k); x(k)  =  tempx + hx; f1  =  feval(fun,x); 
x(k)  =  tempx; J(:,k)  =  (f1-f0)/hx; 
end 
end 