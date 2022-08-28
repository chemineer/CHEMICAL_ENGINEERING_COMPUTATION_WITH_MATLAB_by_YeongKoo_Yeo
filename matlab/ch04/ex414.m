% mgpar.m 
x1 = [0.0932 0.1248 0.1757 0.2000 0.2626 0.3615 0.4750 0.5555 0.6718]; 
Ge = [-0.064 -0.086 -0.120 -0.133 -0.171 -0.212 -0.248 -0.252 -0.245]; 
Y = Ge'; X = [x1.^2.*(1-x1); x1.*(1-x1).^2]'; 
A = inv(X'*X)*X'*Y; A12 = A(2); A21 = A(1);  % Margules parameters 
fprintf('A12 = %g, A21 = %g\n', A12, A21); 
f = @(x) exp((1-x)^2*(A12+2*(A21-A12)*x)) - exp(x^2*(A21+2*(A12-A21)*(1-x))); 
x0 = 0.5; x = fsolve(f,x0)  % find x1 such that gamma1=gamma2