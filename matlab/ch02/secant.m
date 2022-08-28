function x = secant(fun,x1,x2) 
% Implements secant method to find root(s) of an equation 
% input 
% fun: function handle that returns f(x) 
% x1,x2: interval limits containing the root 
% Output
% x: zero of f(x) 
tol = 1e-8; kmax = 1e4; f1 = feval(fun,x1); f2 = feval(fun,x2); 
if f1 == 0, x = x1; return; end 
if f2 == 0, x = x2; return; end 
for k = 1:kmax    
x3 = x2 - f2*(x2 - x1)/(f2 - f1); f3 = feval(fun,x3);    
if abs(f3) <= tol, break; end    
if k >= kmax, disp('The secant method not converged.'); break; end    
x1 = x2; x2 = x3; f1 = f2; f2 = f3;  
end 
x = x3; fprintf('Number of iterations by secant: %g\n', k); 
end 
 