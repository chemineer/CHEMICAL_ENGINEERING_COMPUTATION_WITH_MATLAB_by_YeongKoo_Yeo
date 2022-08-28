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
