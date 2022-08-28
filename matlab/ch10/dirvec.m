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
