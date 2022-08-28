function [xopt,fopt, basic] = barnslp(A,b,c,tol)
% barnslp.m: solution of LP by using Barnes' interior point method
% Problem: Minimize c'x subject to Ax = b (assumed to be non-degenerate)
% Inputs:
% A,b: coefficients of constraint equations
% c: coefficients of the objective function
% Outputs:
% xopt: solution vector
% fopt: function value at xsol
% basic: list of basic variables
% Example: minimize f(x)=2x1 + 3x2 + 2x2 - x4 + x5
% subject to 3x1-3x2+4x3+2x4-x5=0, x1+x2+x3+3x4+x5=2
% A = [3 -3 4 2 -1;1 1 1 3 1]; b = [0; 2]; c = [2 3 2 -1 1]; tol = 1e-6;
% [xopt,fopt, basic] = barnslp(A,b,c,tol)
% Initialization
x2 = [ ]; x = [ ]; [m n] = size(A); ctemp = zeros(1,n); clen = length(c);
for j = 1:clen, ctemp(j) = c(j); end
c = ctemp; aplus1 = b - sum(A(1:m,:)')'; cplus1 = 1000000;
A = [A aplus1]; c = [c cplus1]; B = [ ]; n = n+1; x0 = ones(1,n)'; x = x0;
alpha = 0.0001; lambda = zeros(1,m)'; iter = 0;
% Main step
while abs(c*x - lambda'*b) > tol
x2 = x.*x; D = diag(x); D2 = diag(x2); AD2 = A*D2;
lambda = (AD2*A')\(AD2*c'); dualres = c' - A'*lambda; normres = norm(D*dualres);
for i = 1:n
if dualres(i)>0, ratio(i) = normres/(x(i)*(c(i)-A(:,i)'*lambda));
else, ratio(i) = inf; end
end
R = min(ratio) - alpha; x1 = x - R*D2*dualres/normres; x = x1; basiscount = 0;
B = [ ]; basic = [ ]; cb = [ ];
for k = 1:n
if x(k)>tol, basiscount = basiscount+1; basic = [basic k]; end
end
% Non-degenerate problem
if basiscount == m
for k = basic, B = [B A(:,k)]; cb = [cb c(k)]; end
primalsol = b'/B'; xopt = primalsol; break;
end
iter = iter + 1;
end
xopt = x(basic); fopt = c*x;
end
