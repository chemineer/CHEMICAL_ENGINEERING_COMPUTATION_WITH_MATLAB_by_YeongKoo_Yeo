function x = newtrap(f,df,x0) 
% Implements Newton-Raphson method to find root(s) of an equation 
% input 
% f: function handle that returns f(x) 
% df: function handle that returns df/dx 
% x0: initial guess 

% Output 
% x: zero of f(x) 
tol = 1e-8; kmax = 1e4; f0 = feval(f,x0); df0 = feval(df,x0); 
if f0 == 0, x = x0; return; end 
for k = 1:kmax    
x1 = x0 - f0/df0; f1 = feval(f,x1); df1 = feval(df,x1);    
if abs(f1) <= tol || abs(x1-x0) <= tol, break; end    
if k >= kmax, disp('The Newton-Raphson method not converged.'); break; end    
x0 = x1; f0 = f1; df0 = df1;  
end 
x = x1; fprintf('Number of iterations by Newton-Raphson: %g\n', k); 
end 