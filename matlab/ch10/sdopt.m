function [xopt,fopt, iter] = sdopt(fun, delfun, x0, alpha0, crit, kmax)
% sdopt.m: steepest descent method
% inputs:
% fun: objective function
% delfun: gradient of fun
% x0: starting point
% alpha0: initial step size
% crit: stopping criterion
% kmax: maximum iterations
% outputs:
% xopt: optimal point
% fopt: function value at the optimal point (=f(xopt))
% iter: number of iterations
% Example:
% fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% delfun = @(x) [-400*x(1)*(x(2)-x(1)^2)+2*(x(1)-1); 200*(x(2)-x(1)^2)];
% x0 = [-1.2 1]; crit = 1e-6; alpha0 = 5; kmax = 1e3;
% [xopt,fopt, iter] = sdopt(fun, delfun, x0, alpha0, crit, kmax)

nfmax = 10; % maximum function calls during a line search
ni = 2; % indicate poor interval reductions before sectioning
h = alpha0; beta = 0.9; xz = x0; f = fun(x0); x1 = 0; f1 = f; fs = f; k = 0; nc = 0; fold = f;
while (1)
k = k + 1;
if k > kmax, disp('Number of maximum iterations exceeded.'); break; end
x = xz; delf = delfun(x);
if abs(norm(delf)) <= crit, break; end
d = -delf'; % steepest descent direction vector
d = d / norm(d);  % row vector
[xs, fs] = quadappx(fun, xz, x, d, x1, f1, h, nfmax, ni, crit);
x1 = 0; f1 = fs; xz = xz + xs*d; h = beta * xs; fpr = fs;
if abs(fpr - fold) < crit
nc = nc + 1;
if nc == ni, break; end; else, nc = 0; end
fold = fpr;
end
iter = k; xopt = x; fopt = fs;
end

function [x2, f2] = quadappx(fun, xz, x, d, x1, f1, h, nfmax, ni, crit)
% Quadratic approximation method for line search
tau = (sqrt(5)-1)/2; x2 = x1 + h; f2 = fun(xz + x2*d);
if f2 < f1
while (1)
h = h / tau; x3 = x2 + h; f3 = fun(xz + x3*d);
if f3 > f2, break;
else, f1 = f2; x1 = x2; f2 = f3; x2 = x3; end
end
else
x3 = x2; f3 = f2;
while (1)
x2 = (1 - tau) * x1 + tau * x3; f2 = fun(xz + x2*d);
if f2 <= f1, break; else, x3 = x2; f3 = f2; end
end
end
sf = 0.05;  % 0 < sf < 0.5
if (x1 >= x2 | x2 >= x3), disp('Incorrect interval.'); return;
elseif (f1 <= f2 | f2 >= f3), disp('Not 3-point pattern.'); return;
end
vs = 0; vc = 0; wc = 0; j = 1;
while j <= nfmax
sold = abs(x3-x1); fmold = (f1+f2+f3)/3.;
if vs == 0
A = (x1-x2)*(x1-x3); B=(x2-x1)*(x2-x3); C=(x3-x1)*(x3-x2);
x4 = (f1*(x2+x3)/A + f2*(x1+x3)/B + f3*(x1+x2)/C)/(f1/A+f2/B+f3/C)/2;
else
if x2 <= (x1+x3)/2, x4 = x2 + (1-tau)*(x3-x2);
else, x4 = x3 - (1-tau)*(x2-x1); end
vs = 0;
end
%safeguard against coincident points
dxs = sf*min(abs(x2-x1), abs(x3-x2));
if abs(x4-x1) < dxs, x4 = x1+dxs;
elseif abs(x4-x3) < dxs, x4 = x3-dxs;
elseif abs(x4-x2) < dxs
if x2 > (x1+x3)/2, x4 = x2-dxs;
else, x4 = x2+dxs; end
end
f4 = fun(xz + x4*d);
if (x4 > x2)
if (f4 >= f2), x3 = x4; f3 = f4;
else, x1 = x2; f1 = f2; x2 = x4; f2 = f4; end
else
if (f4 >= f2), x1 = x4; f1 = f4;
else, x3 = x2; f3 = f2; x2 = x4; f2 = f4; end
end
snew = abs(x3-x1); fmnew = (f1+f2+f3)/3.;
if abs(x3-x1) <= crit, break; end
if abs(fmnew-fmold) <= crit
wc = wc + 1;
if wc == 2, break; end; else, wc = 0; end
if snew/sold > tau, vc = vc + 1;
if vc == ni, vc = 0; vs = 1; end
else, vc = 0; vs = 0; end
end
end
