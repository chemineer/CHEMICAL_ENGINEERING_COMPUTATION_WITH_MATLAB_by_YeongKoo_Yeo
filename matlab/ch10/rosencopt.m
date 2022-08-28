function [xopt,fopt,iter] = rosencopt(fun,delfun,x0,A,b,ne,m,crit)
% rosencopt.m: Rosen's gradient projection method
% Problem type:
% Minimize f(x)
% subject to Aj*x = b (j=1,...,ne), Aj*x <= b (j=ne+1,...,m)
% Starting point should satisfy all constraints.
% inputs:
% fun: objective function
% delfun: gradient of fun
% x0: starting point (must satisfy all constraints)
% A: constraints coefficient matrix
% b: right-hand side of constraints (column vector)
% ne: number of equality constraints
% m: number of inequality constraints
% crit: stopping criterion
% outputs:
% xopt: optimal point
% fopt: function value at the optimal point (=f(xopt))
% iter: number of iterations
% Example:
% fun = @(x) -x(1)*x(2)*x(3);
% delfun = @(x) [-x(2)*x(3); -x(1)*x(3); -x(1)*x(2)];
% x0 = 10*[1 1 1]; crit = 1e-4; ne = 0; m = 8;
% A = [-1 0 0;0 -1 0;0 0 -1;1 0 0;0 1 0;0 0 1;-1 -2 -2;1 2 2];
% b = [0 0 0 42 42 42 0 72]';
% [xopt,fopt,iter] = rosencopt(fun,delfun,x0,A,b,ne,m,crit,kmax)
x = x0; n = length(x0); % n: number of variables
f = fun(x); iter = 0;
while (1)
iter = iter + 1; [nc, nca] = actcont(n,m,ne,crit,x,A,b); % active constraints
[dv,d,Bm] = dirvec(delfun,n,nc,ne,crit,x,nca,A); % direction (row) vector
if (abs(dv) < crit), break; end
alphak = maxstep(delfun,n,m,nc,x,x0,d,A,b,nca); % maximum step size
x = x0 + alphak*d; fp = delfun(x)'*d';
if (fp > 0), x = bisec(delfun,n,alphak,x0,d); end
f = fun(x); x0 = x;
end
xopt = x; fopt = f;
end

function [nc, nca] = actcont(n,m,ne,crit,x,A,b)
g = A*x' - b;
for k = 1:m, nca(k) = k; end
nc = ne;
for j = ne+1:m
if (abs(g(j)) < crit), nc = nc + 1; ntemp = nca(j); nca(j) = nca(nc); nca(nc) = ntemp; end
end
end
function x = bisec(delfun,n,alphak,x0,d)
mcrit = 1e-6; a1 = 0; a2 = alphak; aw = a2 - a1;
while (a2 - a1) > mcrit*aw
am = (a1 + a2)/2; x = x0 + am*d; fp = delfun(x)'*d';
if (fp < 0), a1 = am; elseif (fp > 0), a2 = am; else, break; end
end
x = x0 + a1*d;
end

function [dn,d,Bm] = dirvec(delfun,n,nc,ne,crit,x,nca,A)
df = delfun(x)'; % row vector
while (1)
if (nc == 0)
d = -df; dn = norm(d); if (dn < crit), break; end; d = d/dn; Bm(1) = 0; return;
else
for j = 1: nc, ic = nca(j);
for jk = 1: nc
jc = nca(jk); Am(j,jk) = 0.;
for k = 1: n, Am(j,jk) = Am(j,jk) + A(ic,k) * A(jc,k); end
end
end
for j = 1: nc
ic = nca(j); Bm(j) = 0.; for k = 1: n, Bm(j) = Bm(j) - A(ic,k) * df(k); end
end
Am = Am(1:nc,1:nc); Bm = Bm(1:nc)'; Bm = inv(Am)*Bm;
for j = 1: n
d(j) = -df(j); for k = 1: nc, kn = nca(k); d(j) = d(j) - A(kn,j) * Bm(k); end
end
end
dn = norm(d);
if (nc == ne & dn <= crit), break; end
if (dn <= crit)
Bmin = Bm(ne+1); imin = ne + 1;
for j = ne + 1: nc, if (Bmin > Bm(j)), Bmin = Bm(j); imin = j; end; end
if (Bmin >= 0.), break;
else, ntemp = nca(imin); nca(imin) = nca(nc); nca(nc) = ntemp; nc = nc - 1; end
else, d = d/dn; break;
end
end
end

function alphak = maxstep(delfun,n,m,nc,x,x0,d,A,b,nca)
nq = 0;
for k = nc + 1: m
cq = nca(k); c = -b(cq); aq = 0.;
for j = 1:n, c = c + A(cq,j)*x(j); aq = aq + A(cq,j)*d(j); end
if (aq ~= 0)
am = -c/aq;
if (am > 0.), nq = nq + 1;
if (nq == 1), alphak = am;
else, if (alphak > am), alphak = am; end
end
end
end
end
if (nq == 0)
alphak = 1;
while (2)
x = x0 + alphak*d; fp = delfun(x)'*d'; if (fp > 0.), break; end; alphak = 2*alphak;
end
end
end
