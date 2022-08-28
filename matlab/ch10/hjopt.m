function [xopt,fopt,iter] = hjopt(fhj,x0,crit)
% hjopt.m: minimization by the Hooke and Jeeves pattern search method
% Inputs:
% fhj: objective functions
% x0: starting point
% crit: stopping criterion
% Outputs:
% xopt: optimal point
% fopt: objective function value at x=xopt
% iter: number of iterations
% Example:
% x0 = [-3 -1 0 1]; crit = 1e-4;
% f = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
% [xopt,fopt,iter] = hjopt(f,x0,crit)
inist = 0.5; fr = 0.125; x = x0; xb = x; n = length(x0); stsize = inist;
nc = 0; % number of contraction steps
nb = 0; % number of base changes
np = 0; % number of pattern moves
fb = fhj(xb); iter = 1; fold = fb;
while (1)
fk = fb; [fk, x] = search(fhj, n, stsize, fk, x);
if (fb-fk > crit)
while (2)
icv = 0;
for j = 1:n, cp = x(j); x(j) = 2*cp - xb(j); xb(j) = cp; end
fb = fk; fk = fhj(x); [fk, x] = search(fhj, n, stsize, fk, x);
if (fb-fk <= crit)
x = xb;
if (nb > 1), if abs(fk - fold) < crit, icv = 1; break; end; end
nb = nb + 1; break;
end
np = np + 1;
end
if (icv == 1), break; end
else
fold = fk; if (stsize < crit), break; end; stsize = fr*stsize; nc = nc + 1;
end
iter = iter + 1;
end
xopt = x; fopt = fk;
end
function [fk, xk] = search(fhj, n, stsize, fk, x)
for k = 1:n
cpt = x(k); x(k) = cpt + stsize; f = fhj(x);
if (f < fk), fk = f;
else
x(k) = cpt - stsize; f = fhj(x);
if (f < fk ), fk = f; else, x(k) = cpt; end
end
end
xk = x;
end
