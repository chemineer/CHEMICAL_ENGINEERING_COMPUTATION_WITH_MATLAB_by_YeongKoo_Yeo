function [xopt,fopt,iter] = nmopt(fun,x0,crit)
% nmopt.m: minimization by the Nelder and Mead's simplex method
% Inputs:
% fun: objective functions
% x0: starting point
% crit: stopping criterion
% Outputs:
% xopt: optimal point
% fopt: objective function value at x=xopt
% iter: number of iterations
% Example:
% x0 = [-3 -1 0 1]; crit = 1e-4;
% f = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
% [xopt,fopt,iter] = nmopt(f,x0,crit)
rf = 1; % reflection parameter
ef = 2; % expansion parameter
cf = 0.5; % contraction parameter
sf = 0.5; % scale parameter
n = length(x0); beta = 1; n1 = n+1; iter = 0;
x = x0; fv(1) = fun(x);
for j = 1:n, pt(1,j) = x(j); end
% Create simplex and evaluate function values
for k = 2:n1
u = zeros(1,n); u(1,k-1) = 1; pt(k,:) = pt(1,:) + beta*u; x = pt(k,:);
fv(k) = fun(x);
end
while (1)
nloop = 0;
while (2)
nloop = nloop + 1; if (nloop == 2), break; end
iter = iter + 1;
% Find highest, lease, and second highest values
[fh, indh] = max(fv); [fl, indl] = min(fv); fs = fv(indl); inds = indl;
for k = 1:n1
if (k ~= indh), if (fs < fv(k)), fs = fv(k); inds = k; end; end
end
% Find the average of n points xm
for k = 1:n
xm(k) = 0;
for j = 1:n1, if (j ~= indh), xm(k) = xm(k) + pt(j,k); end; end
xm(k) = xm(k)/n;
end
% Reflection procedure
xr = xm + rf*(xm - pt(indh,:)); fr = fun(xr);
if (fr >= fl & fr <= fs) % accept reflection
fv(indh) = fr; pt(indh,:) = xr; break;
end
% Expansion procedure
if (fr < fl )
xe = xr + ef*(xr - xm); fe = fun(xe);
if (fe < fl) % accept expansion
fv(indh) = fe; pt(indh,:) = xe; break;
else % accept reflection
fv(indh) = fr; pt(indh,:) = xr; break;
end
end
% Contraction procedure
if (fr > fh)
xc = xm + cf*(pt(indh,:) - xm); [fc] = fun(xc);
if (fc <= fh) % accept contraction
fv(indh) = fc; pt(indh,:) = xc; break;
end
elseif (fr > fs & fr <= fh)
xc = xm + cf*(xr - xm); fc = fun(xc);
if (fc <= fr) % accept contraction
fv(indh) = fc; pt(indh,:) = xc; break;
end
end
% Scaling
for k = 1:n1
if (k ~= indl)
for j = 1:n, x(j) = sf *pt(k,j) + (1 - sf) *pt(indl,j);
pt(k,j) = x(j); end
fv(k) = fun(x);
end
end
end
sigma = std(fv); avg = mean(fv);
if (sigma <= crit), inc = inc + 1; if (inc == 2), break; end
else, inc = 0; end
end
xopt = pt(indl,:); fopt = fl;
end
