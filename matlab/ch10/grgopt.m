function [xopt,fopt,iter] = grgopt(grgfun,delgrgf,delgrgg,x0,xl,xu,kmax,crit)
% grgopt.m: minimization by the generalized reduced gradient (GRG) method
% Problem type: min. f(x)
% subject to gj(x) = 0 (j=1,...,nca), xl <= x <= xu
% Inputs:
% grgfun: objective function and constraints
% delgrgf: gradient of f
% delgrgg: gradient of g
% x0: starting point (must satisfy all constraints)
% xl,xu: lower and upper limit of x
% kmax: maximum possible iterations
% crit: stopping criterion
% Outputs:
% xopt: optimal point
% fopt: function value at the optimal point (=f(xopt))
% iter: number of iterations
% Example:
% (objective function and constraints are defined by function fg1)
% df1 = @(x) [-1.2 -3 0]; dg1 = @(x) [2*x(1) 12*x(2) 1];
% x0 = [0 0 1]; xl = [0 0 0]; xu = [10 10 10]; crit = 1e-4; kmax = 1e3;
% [xopt,fopt,iter] = grgopt(@fg1,df1,dg1,x0,xl,xu,kmax,crit)

dcrit = crit/10; fcrit = crit*1e-2; x = x0; [f g] = grgfun(x); [ncs nv] = size(delgrgg(x));
for k = 1: ncs, if (abs(g(k)) > crit), disp('Infeasible starting point.');
break; end; end
indc = 0; fold = f;
for iter = 1:kmax
df = delgrgf(x); dg = delgrgg(x); for k = 1:nv, nbr(k) = k; end
[nbr,df,dg] = redgr(nbr,nv,ncs,x,xl,xu,df,dg, crit); % basic variables
dm = 0.;
for k = ncs+1:nv
nk = nbr(k); d(nk) = -df(nk);
if (x(nk) < xl(nk)+crit & df(nk) > 0), d(nk) = 0; end
if (x(nk) > xu(nk)-crit & df(nk) < 0), d(nk) = 0; end
if (dm < abs(d(nk))) dm = abs(d(nk)); end
end
if (dm < crit), break; end
for k = 1: ncs
nk = nbr(k); d(nk) = 0;
for j = ncs+1:nv, nj = nbr(j); d(nk) = d(nk) - dg(k,nj)*d(nj); end
end
x0 = x; alpha = findstep(nv,x,xl,xu,d); x = x0;
for k = 1: nv, xtemp(k) = x(k) + alpha*d(k); end
[ftemp g] = grgfun(xtemp); df = delgrgf(xtemp); dg = delgrgg(xtemp);
fup = 0;
for k = 1:nv, fup = fup + df(k)*d(k); end
if (fup > 0 | ftemp > f)
c = 0; x0 = x; [alpha,ftemp] = goldsec(nv,ncs,c,alpha,dcrit,d,x); x = x0;
for k = 1:nv, xtemp(k) = x(k) + alpha*d(k); end
end
[c g] = grgfun(xtemp); ag = abs(g(1));

for k = 1:ncs, if (ag < abs(g(k))), ag = abs(g(k)); end; end
% If infeasible, apply Newton's method
if (ag > crit)
blim = 20;
for bi0 = 1:blim
bnew = 0; bi2 = 12;
for bi1 = 1:bi2
bnew = bnew + 1; df = delgrgf(xtemp); dg = delgrgg(xtemp);
[g] = rednewton (dg, g, nbr, ncs); bf = 0;
for k = 1: ncs
nk = nbr(k); xtemp(nk) = xtemp(nk) - g(k);
if ((xtemp(nk) < xl(nk)) | (xtemp(nk) > xu(nk))), bf = 1; break; end
end
if (bf == 1), break; end; [c g] = grgfun(xtemp); ag1 = abs(g(1));
for k = 1: ncs
if (ag1 < abs(g(k))), ag1 = abs(g(k)); end
end
if (bf == 0 & ag1 < crit), break; end
bf = 1;
if (bnew > 3 | ag1 > ag), break; end
ag = ag1;
end
if (bf == 0)
[ftemp, g] = grgfun(xtemp); if (ftemp < f), break; end
end
alpha = alpha/2;
for k = 1:nv, xtemp(k) = x(k) + alpha*d(k); end
[c g] = grgfun(xtemp); ag = abs(g(1));
for k = 1: ncs, if (ag < abs(g(k))), ag = abs(g(k)); end; end
end
end
if (bf==0 & ftemp<f)
for k = 1: nv, x(k) = xtemp(k); end; f = ftemp;
end
fcrit = abs(f)*fcrit + fcrit;
if (abs(f - fold) < fcrit)
indc = indc + 1; if (indc > 1), break; end
else, indc = 0; end
fold = f;
end
xopt = x; fopt = f;
end

function [x4,ft] = goldsec(nv,ncs,x1,x4,dcrit,d,x)
tau = (sqrt(5)-1)/2; x2 = tau*x1 + (1-tau)*x4;
for k = 1:nv, xtemp(k) = x(k) + x2*d(k); end
[f2, g] = grgfun(xtemp);
for count = 1:100
x3 = tau*x4 + (1-tau)*x1;
for k = 1:nv, xtemp(k) = x(k) + x3*d(k); end
[f3, g] = grgfun(xtemp);
if (f2 < f3), x4 = x1; x1 = x3;
else, x1 = x2; x2 = x3; f2 = f3; end
if (abs(x4 - x1) <= dcrit), break; end
end

x4 = x2; ft = f2;
end

function [g] = rednewton(dgr, g, nbr, ncs)
for k = 1: ncs-1
nk = nbr(k);
for jk = k+1:ncs
c = dgr(jk,nk)/dgr(k,nk);
for j = k + 1: ncs, nj = nbr(j); dgr(jk,nj) = dgr(jk,nj) - c*dgr(k,nj); end
g(jk) = g(jk) - c*g(k);
end
end
g(ncs) = g(ncs)/dgr(ncs,nbr(ncs));
for jm = 1:ncs-1
jk = ncs - jm; im = nbr(jk); c = 1/dgr(jk,im); g(jk) = c*g(jk);
for k = jk + 1: ncs, nk = nbr(k); g(jk) = g(jk) - c*dgr(jk,nk)*g(k); end
end
end

function [nbr,df,dg] = redgr(nbr,nv,ncs,x,xl,xu,df,dg,xcrit)
for k = 1: ncs
nk = nbr(k); kcount = 0;
for j = k: nv
nj = nbr(j);
if (x(nj) > xl(nj)+xcrit & x(nj) < xu(nj)-xcrit)
kcount = kcount + 1;
if (kcount == 1), pivot = dg(k,nj); jpv = j;
else
if (abs(pivot) < abs(dg(k,nj))), pivot = dg(k,nj); jpv = j; end
end
end
end
nbr(k) = nbr(jpv); nbr(jpv) = nk; nk = nbr(k);
if (k == ncs), break; end
for jk = k+1: ncs
dgratio = dg(jk,nk)/dg(k,nk);
for j = k+1:nv, nj = nbr(j); dg(jk,nj) = dg(jk,nj) - dgratio*dg(k,nj); end
end
end
nbs = nbr(ncs);
for j = ncs + 1: nv
nj = nbr(j); dg(ncs,nj) = dg(ncs,nj)/dg(ncs,nbs);
for jm = 1:ncs-1
jk = ncs-jm; km = nbr(jk); dgratio = 1/dg(jk,km); dg(jk,nj) = dgratio*dg(jk, nj);
for k = jk+1: ncs, nk = nbr(k); dg(jk,nj) = dg(jk,nj) - dgratio*dg(jk,nk)*dg(k,nj); end
end
end
% Calculate reduced gradient
for jk = ncs+1:nv
km = nbr(jk); for j = 1: ncs, nj = nbr(j); df(km) = df(km) - dg(j,km)*df(nj); end
end
end

function [alpha] = findstep(nv, x, xl, xu, d)
% Calculate maximum step size
kcount = 0;
for j = 1:nv
au = xu(j) - x(j); % check upper bound
if (d(j) > 1E-30*(au+1))
kcount = kcount + 1; ad = au/d(j);
if (kcount == 1), alpha = ad; end
if (kcount > 1), if (alpha > ad), alpha = ad; end; end
end
al = x(j) - xl(j); % check lower bound
if (-d(j) > 1E-30*(al + 1))
kcount = kcount + 1; ad = -al/d(j);
if (kcount == 1), alpha = ad; end
if (kcount > 1), if (alpha > ad), alpha = ad; end; end
end
end
end
