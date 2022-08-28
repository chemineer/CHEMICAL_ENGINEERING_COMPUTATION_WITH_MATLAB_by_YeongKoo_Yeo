function [xopt,fopt,iter] = saopt(fun,T,xl,xu,rp,rs,crit)
% saopt.m: minimization by simulated annealing (SA) method
% Inputs:
% fun: objective function
% T: initial temperature
% xl,xu: lower and upper bounds on x
% rp: temperature reduction factor
% rs: step reduction parameter
% crit: stopping criterion
% Outputs:
% xopt: optimal point
% fopt: objective function value at x=xopt
% iter: number of cycle iterations
% Example:
% f = @(x) -cos(5*sqrt(sum(x-5).^2)) + 0.1*sum(x-5).^2;
% T = 100; xl = -10*[1 1]; xu = 10*[1 1]; rp = 0.8; rs = 0.9; crit = 1e-8;
% [xopt,fopt,iter] = saopt(f,T,xl,xu,rp,rs,crit)
% Initialization
sf = 2; stp = 1; np = 10; nc = 20; nt = 2e4; ir = 16; rcrit = 1e-10; n = length(xu);
% Feasible starting point
for j = 1:n, x(j) = xl(j) + rand*(xu(j) - xl(j)); xs(j) = x(j); xmin(j) = x(j); end
f = fun(x); fmin = f; fold = f; citer = 0;
% Set step sizes, step factors and acceptance ratios
for j = 1:n, stsize(j) = stp; ar(j) = 1; end
while (1)
for piter = 1:np % temperature step loop
for ic = 1:nc % % search cycles: search along coordinate direction
for k = 1:n
xs(k) = x(k) + (2*rand - 1)*stsize(k);
if (xs(k) < xl(k)) | (xs(k) > xu(k)), xs(k) = xl(k) + rand*(xu(k) - xl(k));
end
fs = fun(xs);
if fs <= f % point is accepted: update xmin and fmin
x(k) = xs(k); f = fs;
if fs < fmin, xmin = xs; fmin = fs; end
else
p = exp((f - fs)/T);
if rand < p, x(k) = xs(k); f = fs;
else % point rejected
xs(k) = x(k); ar(k) = ar(k) - 1/nc;
end
end
end
end
% Adjust step so that about half the points are accepted.
for j = 1:n
if ar(j) > 0.6, stsize(j) = stsize(j)*(1 + sf*(ar(j) - 0.6)/0.4);
elseif ar(j) < 0.4, stsize(j) = stsize(j)/(1 + sf*(0.4 - ar(j))/0.4); end
if stsize(j) > xu(j) - xl(j), stsize(j) = xu(j) - xl(j); end
ar(j) = 1;
end
end
stsize = stp*(xu - xl); fcrit = crit + rcrit*abs(fmin);
if (fmin <= fold) & (fold - fmin < fcrit)
citer = citer + 1; if citer >= ir, break; end
else, citer = 0;
end
% Reduce temperature
T = rp*T; stp = rs*stp; x = xmin; f = fmin; fold = f;
end
xopt = xmin; fopt = fmin; iter = citer;
end
