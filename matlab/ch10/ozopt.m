function [xopt,fopt,iter] = ozopt(A,c,nl,ne)
% ozopt.m: zero-one programming method for inter minimization problem
% inequality constraints should be rearranged to <= constraints.
% Inputs:
% A: coefficient matrix for constraints including right hand sides
% c: coefficient vector for objective function (row vector)
% nl: number of <= constraints
% ne: number of = constraints
% Outputs:
% xopt: optimal point (0 or 1)
% fopt: function value at the optimal point (=f(xopt))
% iter: number of iterations
% Example
% nl = 3; ne = 0;
% A = [-1 3 6 6; -2 -3 -3 -2; -1 0 -1 -1]; c = [4 5 3];
% [xopt,fopt,iter] = ozopt(A,c,nl,ng,ne)
% nl = 3; ne = 3;
% A = [1 0 0 1 0 1 0 0 0 1;0 1 0 0 0 0 1 0 0 1;0 0 0 1 0 0 0 0 1 1;1 1 1 0 0 0 0 0 0 1;0 0 0 1 1 0 0 0 0 1;0 0 0 0 0 0 0 1 1 1];
% c = -[5 3 1 3 5 2 5 5 2];
nv = length(c); m = nl + ne; x = zeros(1,nv); Am = zeros(m+1,nv+1); xvec = zeros(1,nv);
indv = zeros(1,nv); temA = zeros(1,m+1); xmin = zeros(1,nv);
[rAm,cAm] = size(A); Am(1:rAm,1:cAm) = A;
for k = 1:m, Am(k,nv+1) = -Am(k,nv + 1); end
% Convert <= constraints to >=
for k = 1:nl, for j = 1:nv+1, Am(k,j) = -Am(k,j); end; end
% Convert = constraints to >= by adding one >= constraint
if ne > 0
m = m + 1;
for k = nl+1:m-1
for j = 1:nv+1, Am(m,j) = Am(m,j) - Am(k,j); end
end
end
% Convert variables if c(j) < 0
check0 = 0;
for j = 1:nv
if c(j) < 0
indv(j) = 1; check0 = check0 + c(j); c(j) = -c(j);
for k = 1:m, Am(k,nv+1) = Am(k,nv+1) + Am(k,j); Am(k,j) = -Am(k,j); end
end
end
for k = 1:nv, xvec(k) = k; x(k) = 0; end
for k = 1:m, temA(k) = Am(k,nv+1); end
indf = 1; inds = 0; f = 0; iter = 0;
while (1)
sflag = 0; iter = iter + 1;
if xvec(indf) > -1
sflag = 1;
for k = 1:m, if temA(k) < 0, sflag = 0; break; end; end
end
if sflag == 1 % feasible solution
inds = inds + 1;
if inds == 1 || (inds > 1 && f < fmin)
fmin = f;
for k = 1:nv, xmin(k) = x(k); end
end
end
cflag = 0; nflag = 0;
if sflag == 0 % infeasible solution
for k = 1:m
indk = temA(k);
if indk < 0
for j = indf+1:nv
if Am(k,xvec(j)) > 0, indk = indk + Am(k,xvec(j)); end
end
end
if indk < 0, cflag = 1; break; end
end
if inds > 0
nflag = 1;
for k = indf+1:nv, if f + c(xvec(k)) < fmin, nflag = 0; break; end; end
end
end
if sflag == 1 || cflag == 1 || nflag == 1
while xvec(indf) < 0
xvec(indf) = -xvec(indf); indf = indf - 1;
if indf == 0, break; end
end
if indf == 0, break; end
x(xvec(indf)) = 0;
% Update constraints and function value (f)
for k = 1:m, temA(k) = temA(k) - Am(k,xvec(indf)); end
f = f - c(xvec(indf)); xvec(indf) = -xvec(indf);
else
indf = indf + 1;
for k = indf:nv
difa = 0;
for j = 1:m
temg = temA(j) + Am(j,xvec(k)); if temg < 0,
difa = difa - temg; end
end
if k == indf, delf = difa; ink = k;
else, if delf > difa, delf = difa; ink = k; end; end
end
inj = xvec(indf); xvec(indf) = xvec(ink); xvec(ink) = inj;
% Update x and f
x(xvec(indf)) = 1;
for k = 1:m, temA(k) = temA(k) + Am(k,xvec(indf)); end
f = f + c(xvec(indf));
end
end
if inds == 0, display('No feasible solution.'); return; end
% Assign outputs
for k = 1:nv, if (indv(k) == 1), xmin(k) = 1 - xmin(k); end; end
xopt = xmin; fopt = fmin;
end
