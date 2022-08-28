function x = bisectn(fun,x1,x2) 
% input 
% fun: function handle that returns f(x) 
% x1,x2: interval limits containing the root 
% Output 
% x: zero of f(x) 
tol = 1e-8; kmax = ceil(log(abs(x2-x1)/tol)/log(2)); 
f1 = feval(fun,x1); f2 = feval(fun,x2); 
if f1 == 0, x = x1; return; end 
if f2 == 0, x = x2; return; end 
if f1*f2 > 0, error('The root is not located in [x1, x2].'); end 
for k = 1:kmax    
x3 = (x1+x2)/2; f3 = feval(fun,x3); 
if abs(f3) <= tol, x = x3; break; end    
if f2*f3 < 0, x1 = x3; f1 = f3; else x2 = x3; f2 = f3; end 
end 
x = (x1+x2)/2; 
fprintf('Number of iterations by bisect: %g\n', k); 
end 