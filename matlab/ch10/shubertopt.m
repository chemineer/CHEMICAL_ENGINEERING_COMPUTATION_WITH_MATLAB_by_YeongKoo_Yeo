function [xopt, fopt, nf] = shubertopt(objfun, a, b, C, crit, nfmax)
% shubertopt.m: optimization by Shubert-Piyavskii algorithm
% 1-D maximization of Lipschitz functions

% inputs:
% objfun: objective function
% a,b: end points of initial interval
% C: Lipschitz constant
% crit: tolerance
% nfmax: maximum number of function evaluations
% outputs:
% xopt: optimum point
% fopt: function value at the optimum point
% nf: number of function evaluations
% Sample run:
% C = 8; a = -3; b = 8; crit = 1e-6; nfmax = 2000;
% fun = @(x) -sin(x)-sin(3.5*x);
% [xopt,fopt,nf] = shubertopt(fun, a, b, C, crit, nfmax)

nf = 0; x0 = (a + b)/2; nf = nf + 1; y0 = objfun(x0); ymax = y0; xmax = x0; fmax = y0 + C*(b - a)/2;
T(1) = b; Z(1) = y0 + C*(b - a)/2; T(2) = a; Z(2) = y0 + C*(b - a)/2; n = 2;
while ((fmax - ymax) > crit & nf <= nfmax)
tn = T(n); zn = Z(n); nf = nf + 1; yn = objfun(tn);
if (yn > ymax), ymax = yn; xmax = tn; end
zL = (zn + yn)/2; zR = zL; tL = tn - (zn - yn)/2/C; tR = tn + (zn - yn)/2/C;
% Replace T(n) and Z(n) by tl,zl and tr,zr
ind1 = 0; ind2 = 0;
if (tL >= a & tL <= b), ind1 = 1; end
if (tR >= a & tR <= b), ind2 = 1; end
if (ind1 == 1 & ind2 == 0), T(n) = tL; Z(n) = zL;
elseif (ind1 == 0 & ind2 == 1), T(n) = tR; Z(n) = zR;
elseif (ind1 == 1 & ind2 == 1), T(n) = tL; Z(n) = zL; T(n+1) = tR; Z(n+1) = zR;
n = n+1; end
[Z, Indx] = sort(Z); T = T(Indx); fmax = Z(n);
end
xopt = xmax; fopt = ymax;
end
