function [xopt,fopt,istag] = rbopt(fros,x0,crit)
% rbopt.m: minimization by the Rosenbrock’s method
% Inputs:
% fros: objective functions
% x0: starting point
% crit: stopping criterion
% Outputs:
% xopt: optimal point
% fopt: objective function value at x=xopt
% istag: number of stages
% Example:
% x0 = [-3 -1 0 1]; crit = 1e-4;
% f = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
% [xopt,fopt,istag] = rbopt(f,x0,crit)
n = length(x0); % number of variables
mi = 7; stsize = 1; alpha = 3; beta = 0.5; xt = x0;
% Initial vectors along coordinates
for j = 1:n
for k = 1:n, v(j,k) = 0; end; v(j,j) = 1;
end
ft = fros(xt); istag = 0;
while (1)
istag = istag + 1;
% Initialization
for j = 1:n, stp(j) = stsize; stj(j) = 0; sfv(j,1) = 0; sfv(j,2) = 0; end
iter = 0; icont = 0; idn = 0;
while (2)
iter = iter + 1; if iter > n, iter = 1; end
for j = 1:n, x(j) = xt(j) + stp(iter)*v(j,iter); end; f = fros(x);
if f < ft
stj(iter) = stj(iter) + stp(iter); % successful
for j = 1:n, xt(j) = x(j); end
ft = f; idn = 0; stp(iter) = alpha*stp(iter);
if sfv(iter,1) == 0, sfv(iter,1) = 1; icont = icont + 1; end
else % failure
stp(iter) = -beta*stp(iter);
if sfv(iter,2) == 0, sfv(iter,2) = 1; icont = icont + 1; end
idn = idn + 1;
end
if icont == 2*n, break; end
if idn > 50, mi = mi + 1; return; end
end
if istag == 1, fold = ft;
else
if abs(ft - fold) < crit, mi = mi + 1; display('Convergence achieved.');
break; end
end
fold = ft;
% Construct new vectors
for k = 1:n
for j = 1:n, u(k,j) = 0; end
end
for i = 1:n
for j = 1:n, for k = 1:n, u(k,i) = u(k,i) + stj(j)*v(k,j); end; end
end
% Orthogonalization by Gram-Schmidt procedure
for i = 1:n
if i > 1
for j = 1:n, v(j,i) = 0; end
for k = 1:i-1
c = 0; for j = 1:n, c = c + u(j,i)*v(j,k); end; for j = 1:n, v(j,i) = v(j,i) + c*v(j,k); end
end
for j = 1:n, u(j,i) = u(j,i) - v(j,i); end
end
c = 0; for j = 1:n, c = c + u(j,i)*u(j,i); end
c = sqrt(c);
% Take new step length as the length of the first vector
if i == 1
stsize = c;
if stsize < crit, mi = mi + 1; break; end
if c < crit % Orthogonality is lost: reset vectors.
for j = 1:n
for k = 1:n, v(j,k) = 0; end
v(j,j) = 1; break;
end
end
for j = 1:n, v(j,i) = u(j,i)/c; end
end
end
end
xopt = xt; fopt = ft;
end
