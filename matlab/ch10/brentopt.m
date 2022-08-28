function [xopt, fopt, nf] = brentopt(objfun, x1, h, crit)
% brentopt.m: minimization(1-dimensional search) by Brent's algorithm
% input:
% objfun: objective function
% x1: initial point
% h: step size
% crit: convergence criterion
% output:
% x: optimum point
% f: function value at the optimum point
% hf: number of function evaluations for search
% Sample run:
% x1 = 0.01; h = 0.2; crit = 1e-6;
% fun = @(x) exp(x) -3*x + 0.02/x - 0.00004/x;
% [xopt, fopt, nf] = brentopt(fun, x1, h, crit)

% Initialization
clarge = 1e40; tau = (sqrt(5) - 1)/2; nf = 1; f1 = objfun(x1);
%  Determine initial 3-points
x2 = x1 + h; nf = nf+1; f2 = objfun(x2);
if f2 > f1, temp = x1; x1 = x2; x2 = temp; temp = f1; f1 = f2; f2 = temp; h = -h; end
while x1 < clarge
h = h/tau; x4 = x2 + h; nf = nf+1; f4 = objfun(x4);
if f4 > f2, break; end; f1 = f2; x1 = x2; f2 = f4; x2 = x4;
end
x3 = x4; f3 = f4; a = x1; b = x3; fa = f1; fb = f3; ha = a; hb = b;
% Check whether ha < hb
if (b < a), ha = b; hb = a; end
x = x2; fx = fb; w = x; v = w; ev = 0; fx = f2; fw = fx; fv = fw;
while 1
hm = (ha + hb)/2;
%  Check interval convergence
if (abs(x - hm) <= crit - (hb - ha)/2), x2 = x; f2 = fx; return; end
etemp = 0; p = 0; q = 0;
if (abs(ev) > crit)
r = (x - w)*(fx - fv); q = (x - v)*(fx - fw); p = (x - v)*q - (x - w)*r; q = 2*(q - r);
if (q > 0), p = -p; else, q = -q; end;
etemp = ev; ev = d; % length of the larger interva
end
ind1 = 1; if ((q * (ha - x) - p) < 0), ind1 = -1; end
ind2 = 1; if ((q * (hb - x) - p) < 0), ind2 = -1; end
if ((abs(p) >= abs(q*etemp/2)) | (ind1 == ind2))
if (x < hm), ev = hb - x; else, ev = ha - x; end; d = (1 - tau)*ev;
else
d = p/q; u = x + d;
if (((u - ha) < crit) | ((hb - u) < crit)), if (x < hm), d = crit; else, d = -crit;
end; end
end
if (abs(d) >= crit), u = x + d;
else, if (d > 0), u = x + crit; else, u = x - crit; end; end
nf = nf+1; fu = objfun(u);
% Set a, b, x, u, v, w for the next iteration
if (fu <= fx)
if (u < x), hb = x; fb = fx;
else, ha = x; fa = fx; end
v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
else
if (u < x), ha = u; fa = fu;
else, hb = u; fb = fu; end
if ((fu <= fw) | (w == x)), v = w; fv = fw; w = u; fw = fu;
elseif ((fu <= fv) | (v == x) | (v == w)), v = u; fv = fu; end
end
% Check function convergence
if ((fa - fx) + (fb - fx) < crit), x2 = x; f2 = fx; return; end
xopt = x2; fopt = f2;
end
end
