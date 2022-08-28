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
