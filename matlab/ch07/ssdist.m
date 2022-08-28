function f = ssdist(x,dpar)
alpha = dpar.alpha; n = dpar.n; nf = dpar.nf; F = dpar.F;
zf = dpar.zf; q = dpar.q; R = dpar.R; D = dpar.D;
Lr = R; B = F - D; Ls = R + F*q; Vs = Ls - B; Vr = Vs + F*(1-q);
y = alpha*x./(1 + (alpha-1)*x);
f(1) = (Vr*y(2) - (D + R)*x(1)); % condenser
for i = 2:nf-1, f(i) = Lr*x(i-1) + Vr*y(i+1) - Lr*x(i) - Vr*y(i); end
f(nf) = Lr*x(nf-1) + Vs*y(nf+1) - Ls*x(nf) - Vr*y(nf) + F*zf;
for i = nf+1:n-1, f(i) = Ls*x(i-1) + Vs*y(i+1) - Ls*x(i) - Vs*y(i); end
f(n) = (Ls*x(n-1) - B*x(n) - Vs*y(n)); % reboiler
end
