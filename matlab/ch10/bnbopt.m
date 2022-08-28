function [xopt,fopt,iter] = bnbopt(A,b,c,nl,ng,ne,ibd)
% bnbopt.m: branch and bound method for mixed inter minimization problem
% Inputs:
% A: coefficient matrix for constraints
% b: right hand side of constraints (row vector)
% c: coefficient vector for objective function (row vector)
% nl: number of <= constraints
% ng: number of >= constraints
% ne: number of = constraints
% ibd: variable index vector
% Outputs:
% xopt: optimal point (0 or 1)
% fopt: function value at the optimal point (=f(xopt))
% iter: number of iterations
critd = 1e-4; critf = 5e-2; kmax = 1e4; ksol = 0; kf = 0; fmin = 1e30;
nv = length(c); nbd = length(ibd);
for j = 1:length(b)
if (b(j) < 0), disp('b(j) must be non-negative.'); return; end
end
for k = 1:nv, vtype(k) = 1; end
for k = 1:nbd, vtype(ibd(k)) = 2; end
for j = 1:nv, kvar(j) = j; end
iter = 0; kconv = 0;
while (1)
iter = iter + 1;
if (iter > kmax), disp('Maximum possible iterations exceeded'); break; end
ktrac = 0; [f,x,kflag] = linpm(nv,nl,ng,ne,kf,kvar,A,b,c);
if (nbd == 0)
if (kflag > 0), disp('No feasible solution.'); return;
else, fmin = f; xmin = x; break; end
end
jc = 0;
while (2)
jc = jc+1; if (jc == 2), break; end
if (kflag > 0)
if (kf == 0), disp('No feasible solution.'); return; end
ktrac = 1; break; % backtrack operation
else
if (iter == 1), flb = f; end % best low bound
if (f >= fmin), ktrac = 1; break; end
vmax = 0;
for j = 1:nv
if(vtype(j) == 2)
diffx = abs(x(j)-round(x(j)));
if (diffx > vmax), vmax = diffx; jv = j; end
end
end
if (vmax < critd)
ksol = 1;
if (f < fmin), fmin = f; for j = 1:nv, xmin(j) = x(j); end; end
gapf = abs(flb - fmin)/abs(flb);
if (gapf <= critf), kconv = 1; break; end
ktrac = 1; % backtrack operation
break
else
kf = kf + 1; jt = kvar(kf); kvar(kf) = kvar(jv); kvar(jv) = jt;
end
end
end
while (ktrac == 1)
if (kvar(kf) > 0), kvar(kf) = -kvar(kf); break;
else
kvar(kf) = -kvar(kf); kf = kf - 1;
if (kf == 0)
if(ksol == 0), disp('No feasible solution.'); return; end
kconv = 1; break;
end
end
end
if (kconv == 1), break; end
end
xopt = xmin; fopt = fmin;
end
function [f,x,kflag] = linpm(nv,nl,ng,ne,kf,kvar,A,B,C)
if (kf > 0)
idn = [1:nv];
for k = 1:kf, jk = abs(kvar(k)); idn(jk) = 0; end
difnv = nv - kf; nej = ne; Aj = A; Bj = B; f0 = 0; cj = 0; jm = 0;
for k = 1:nv
if idn(k) > 0, cj = cj+1; Ck(cj) = C(k);
else
jm = jm+1;
if (kvar(jm) > 0), cm = 1; else, cm = 0;
end
f0 = f0 + C(k)*cm;
end
end
for k = 1:nl, conk(k) = -1; end
for k = 1:ng, conk(nl+k) = 1; end
for k = 1:nl+ng+ne
cj = 0; jm = 0;
for j = 1:nv
if (idn(j)) > 0, cj = cj + 1; Aj(k,cj) = A(k,j);
else
jm = jm+1;
if (kvar(jm) > 0), cm = 1; else, cm = 0; end
Bj(k) = Bj(k) - A(k,j)*cm;
end
end
if (k <= nl+ng & Bj(k) < 0), conk(k) = -conk(k); Aj(k,:) = -Aj(k,:); Bj(k) = -Bj
(k); end
end
nlk = 0; ngk = 0;
for k = 1:nl+ng
if (conk(k) == -1), nlk = nlk + 1;
elseif (conk(k)==1), ngk = ngk + 1;
end
end
[ids,xk] = sort(conk); A = Aj; B = Bj;
for k = 1:nl+ng, A(k,:) = Aj(xk(k),:); B(k) = Bj(xk(k)); end
nk = nlk + ngk + nej; Aj = A(1:nk,1:difnv); Bj = B(1:nk);
else
difnv = nv; nlk = nl; ngk = ng; nej = ne; Aj = A; Bj = B; Ck = C;
end
[fk,xs,kflag] = linsimx(difnv,nlk,ngk,nej,Aj,Bj,Ck);
if (kflag > 0), x = xs; f = fk; return; end
if (kf > 0)
cj = 0; jm = 0;
for j = 1:nv
if (idn(j)) > 0. cj = cj+1; x(j) = xs(cj);
else
jm = jm+1;
if (kvar(jm)>0), cm = 1; else, cm = 0; end
x(j) = cm;
end
end
f = fk + f0;
else
f = fk; x = xs; return;
end
end

function [f,x,kflag] = linsimx(nv,nl,ng,ne,A,B,C)
% linear programming by simplex method
mnp = 20; bigm = 0; kflag = 0;
nc = nv + nl + ne + 2*ng; nr = nl + ne + ng + 1;
% Initialization
Aj(1:nr,1:nc) = 0; Bj(1:nr) = 0; Aj(1:nr-1,1:nv) = A(1:nr-1,1:nv);
Aj(nr,1:nv) = C(1:nv); Bj(1:nr-1) = B(1:nr-1); A = Aj; B = Bj;
for j = 1:nv, A(nr,j) = C(j); end
for j = 1:nv, bigm = bigm + mnp*abs(A(nr,j)); end
if (nl > 0) % slack variables
for k = 1:nl, A(k,nv+k) = 1; Bs(k) = nv + k; end
end
if (ng > 0)
for k = 1:ng
A(nl+k,nv+nl+k) = -1; A(nl+k, nv+nl+ng+k) = 1; Bs(nl+k) = nv + nl + ng + k;
end
end
if (ne > 0)
for k = 1:ne, A(nl+ng+k,nv+nl+2*ng+k) = 1; Bs(nl+ng+k) = nv + nl + 2*ng + k; end
end
mge = ng + ne;
if (mge > 0), for k = 1:mge, A(nr,nv+nl+ng+k) = bigm; end
end
if (mge > 0) % Remove artificial variables
for k = 1:mge
C = A(nr,nv+nl+ng+k); for j = 1:nc, A(nr,j) = A(nr,j) - C*A(nl+k,j); end
B(nr) = B(nr) - C*B(nl+k);
end
end

% Simplex method
while (1)
jflag = 0;
for j = 1:nc, if (A(nr,j) < 0), jflag = 1; break; end; end
if (jflag == 0), break; end
C = bigm;
for j = 1:nc, if (A(nr,j) < C), C = A(nr,j); idv = j; end; end
jn = 0; jk = 0;
for k = 1:nr-1
if (A(k,idv) > 0)
jk = jk + 1; Dm = B(k)/(A(k,idv) + 1e-10);
if (jk == 1), C = Dm; jp = k;
else, if (Dm < C), C = Dm; jp = k; end
end
jn = 1;
end
end
Bs(jp) = idv; if (jn == 0), kflag = 1; f = 0; x = 0; return; end
Dm = 1/A(jp,idv); B(jp) = Dm*B(jp);
for j = 1:nc, A(jp,j) = Dm*A(jp,j); end
for k = 1:nr
if (k ~= jp)
Em = A(k,idv); for j = 1:nc, A(k,j) = A(k,j) - Em*A(jp,j); end
B(k) = B(k) - Em*B(jp);
end
end
end
x(1:nv) = 0; nt = nv + nl + ng;
for k = 1:nr-1, if (Bs(k) > nt), kflag = 1; f = 0; x = 0; return; end; end
for k = 1:nv
for j = 1:nr-1, if (Bs(j) == k), x(k) = B(j); break; end; end
end
f = B(nr); f = -f; % minimization
end
