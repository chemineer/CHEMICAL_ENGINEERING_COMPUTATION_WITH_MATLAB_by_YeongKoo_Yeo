function [xopt,fopt,iter] = sqpopt(fun,dfun,x0,lam0,mu0,crit)
% sqpopt.m: minimization by using SQP algorithm
% Problem type: minimize f(x) subject to h(x) = 0, g(x) >= b
% Inputs:
% fun: objective and constraint functions
% dfun: gradients of the objective and constraint functions
% x0: starting point
% lam0,mu0: initial Lagrange multipliers
% crit: stopping criterion
% Outputs:
% xopt: optimal point
% fopt: objective function value at x=xopt
% iter: number of iterations
% Example:
% x0 = [3 3]'; lam0 = 1; mu0 = 4; crit = 1e-6;
% [xopt,fopt,iter] = sqpopt(@fun,@dfun,x0,lam0,mu0,crit)
x = x0(:); n = length(x); % number of variables
lam1 = lam0 + 1; Hj = eye(n); fv = fun(x); aj = fv(2:lam1); cj = fv((lam1+1):(lam1+mu0));
Gj = dfun(x); gj = Gj(:,1); Aej = Gj(:,2:lam1)'; Aij = Gj(:,(lam1+1):(lam1+mu0))'; iter = 0; d = 1;
while d >= crit
delx = quadpr(Hj,gj,Aej,-aj,Aij,-cj,zeros(n,1),crit);
ad = Aij*(x+delx) + cj; k = find(ad <= crit); nk = length(k); muj = zeros(mu0,1);
if nk == 0, lamj = inv(Aej*Aej')*Aej*(Hj*delx+gj);
else
Aaik = Aij(k,:); Aaj = [Aej; Aaik]; mun = inv(Aaj*Aaj')*Aaj*(Hj*delx+gj);
lamj = mun(1:lam0); mujh = mun(lam1:end); muj(k) = mujh;
end
alpha = linsearch(fun,x,delx,lam1,muj,crit); delx = alpha*delx; x = x + delx;
grd = dfun(x);
grd1 = grd(:,1); Agrd = grd(:,2:lam1)'; Am = grd(:,(lam1+1):(lam1+mu0))';
gamj = (grd1-gj)-(Agrd-Aej)'*lamj-(Am-Aij)'*muj; qj = Hj*delx;
dg = delx'*gamj; dq = delx'*qj;
if dg >= 0.2*dq, theta = 1;
else, theta = 0.8*dq/(dq-dg); end
eta = theta*gamj + (1-theta)*qj; invdq = 1/dq; invde = 1/(delx'*eta);
Hj = Hj + invde*(eta*eta') - invdq*(qj*qj'); Aej = Agrd; Aij = Am; gj = grd1;
fv = fun(x);
aj = fv(2:lam1); cj = fv((lam1+1):(lam1+mu0)); d = norm(delx); iter = iter + 1;
end
xopt = x; fopt = fv(1);
end

function alpha = linsearch(fun,xj,dx,lam,muj,crit)
% Determine step size by line search method
nmuj = length(muj); alrange = 0:0.01:1; nint = length(alrange);
hz = zeros(nint,1);
for j = 1:nint
xdj = xj + alrange(j)*dx; fv = fun(xdj); af = fv(2:lam);
cf = fv((lam+1):(lam+nmuj));
hz(j) = fv(1) + 1e2*sum(af.^2)- muj'*cf;
end
[mval,indhz] = min(hz); atemp = alrange(indhz); indmu = find(muj <= crit); mc = length(indmu);
if mc == 0, alpha = 0.95*atemp;
else
dv = zeros(mc,1);
for k = 1:mc
for j = 1:nint
aj = alrange(j); xdj = xj + aj*dx; fcj = fun(xdj);
cj = fcj((lam+1):(lam+nmuj)); hz(j) = cj(indmu(k));
end
indhz = find(hz < 0); hc = length(indhz);
if hc == 0, dv(k) = 1;
else, dv(k) = alrange(indhz(1)-1); end
end
mdv = min(dv); alpha = 0.95*min(atemp,mdv);
end
end

function xqr = quadpr(Q,c,Aeq,beq,Ane,bne,x0,crit)
% minimization by quadratic programming method
nbeq = length(beq); nbne = length(bne);
rnbne = nbne + 1.5*sqrt(nbne); aone = 1-crit; x = x0(:); y = Ane*x - bne;
zeta = zeros(nbeq,1); gama = ones(nbne,1); ym = y.*gama; sumy = sum(ym); iter = 0;
while sumy > crit
tau = sumy/rnbne; resid = -Q*x - c + Aeq'*zeta + Ane'*gama;
diffr = beq - Aeq*x; numt = tau - y.*gama; yj = gama./y; ysj = numt./y;
ydm = diag(yj);
Gr = inv(Q + Ane'*ydm*Ane); ag = Aeq*Gr*Aeq';
ayj = resid + Ane'*ysj; diffag = diffr - Aeq*Gr*ayj;
delz = inv(ag)*diffag;
dx = Gr*(ayj + Aeq'*delz); dy = Ane*dx; dgama = (numt - (gama.*dy))./y;
indny = find(dy < 0); yr = min(y(indny)./(-dy(indny)));
indny = find(dgama < 0); gamr = min(gama(indny)./(-dgama(indny)));
aj = aone*min([1 yr gamr]); x = x + aj*dx; gama = gama + aj*dgama;
zeta = zeta + aj*delz;
y = Ane*x - bne; ym = y.*gama; sumy = sum(ym);
iter = iter + 1;
end
xqr = x;
end
