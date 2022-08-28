function dx = dyndist(t,x,dpar,dels)
% System of differential equations for a binary distillation column
% Equilibrium relation: y=alpha*x/(1+(alpha-1)*x)
% Solver: [t x] = ode45(@dyndist,[t0 tf],x0,[],dpar,dels)
% Parameters
alpha = dpar.alpha; n = dpar.n; nf = dpar.nf; Fi = dpar.F;
zfi = dpar.zf; q = dpar.q; Ri = dpar.R; Vsi = dpar.Vs;
md = dpar.md; mb = dpar.mb; mt = dpar.mt;
delR = dels.delR; delRt = dels.delRt; delV = dels.delV; delVt = dels.delVt;
delz = dels.delz; delzt = dels.delzt; delF = dels.delF; delFt = dels.delFt;
% Changes in operating conditions
if t < delRt, R = Ri; else R = Ri + delR; end
if t < delVt, Vs = Vsi; else Vs = Vsi + delV; end
if t < delzt, zf = zfi; else zf = zfi + delz; end
if t < delFt, F = Fi; else F = Fi + delF; end
% Floe rates
Lr = R; Ls = R + F*q; B = Ls - Vs; D = F - B; Vr = Vs + F*(1-q);
% Initialization and phase equilibrium
dx = zeros(n,1); y = alpha*x./(1 + (alpha-1)*x);
% Condenser
dx(1) = Vr*(y(2)-x(1))/md;
% Rectifying section
for i = 2:nf-1, dx(i) = (Lr*x(i-1)+Vr*y(i+1)-Lr*x(i)-Vr*y(i))/mt; end
% Feed stage
dx(nf) = (Lr*x(nf-1)+Vs*y(nf+1)+F*zf-Ls*x(nf)-Vr*y(nf))/mt;
% Stripping section
for i = nf+1:n-1, dx(i) = (Ls*x(i-1)+Vs*y(i+1)-Ls*x(i)-Vs*y(i))/mt; end
% Reboiler
dx(n) = (Ls*x(n-1)-B*x(n)-Vs*y(n))/mb;
end
