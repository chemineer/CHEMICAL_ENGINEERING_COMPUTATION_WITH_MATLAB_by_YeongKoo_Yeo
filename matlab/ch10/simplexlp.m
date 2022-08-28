function [xopt, fopt] = simplexlp(A, b, c, constr)
% simplexlp.m: 2-phase LP minimization problem
% Problem: Minimize f(x) = c*x subject to Ax constr b, x >= 0
% inputs:
% A,b: coefficients of constraint equations (and inequalities)
% c: coefficients of the objective function
% constr: type of constraints (example: constr = '==<>>=')
% outputs:
% xopt: optimum point
% fopt: value of the objective function at x=xopt
% Example: Minimize f(x)=2x1 + 3x2 + 2x2 - x4 + x5
% subject to 3x1-3x2+4x3+2x4-x5=0, x1+x2+x3+3x4+x5=2
% A = [3 -3 4 2 -1;1 1 1 3 1]; b = [0; 2]; c = [2 3 2 -1 1]; constr = '==';
% [xopt, fopt] = simplexlp(A, b, c, constr)
b = b(:); c = c(:)'; [m, n] = size(A); n1 = n; nleq = 0; neq = 0; ncomv = 0;
if length(c) < n, c = [c zeros(1,n-length(c))]; end
for j = 1:m
temx = zeros(m,1); temx(j) = 1;
if(constr(j) == '<') % <=: less than or equal to
A = [A temx]; nleq = nleq + 1;
elseif(constr(j) == '>') % >=: greater than or equal to
A = [A -temx];
else % =: equality constraints
neq = neq + 1;
end
end
lenA = length(A);
if nleq == m
c = [c zeros(1,lenA-length(c))]; A = [A;c]; A = [A [b;0]];
[maux, A, z] = compsim(A, n1+1:lenA, 1, 1);
else
A = [A eye(m) b];
if m > 1, w = -sum(A(1:m,1:lenA)); else, w = -A(1,1:lenA); end
c = [c zeros(1,length(A)-length(c))]; A = [A;c]; 
A = [A;[w zeros(1,m) - sum(b)]];
maux = lenA+1:lenA+m; mv = maux; [maux, A, z] = compsim(A, maux, 2, 1);
nc = lenA + m + 1; x = zeros(nc-1, 1); x(maux) = A(1:m, nc); xm = x(mv);
incomv = intersect(maux, mv);
if (any(xm) ~= 0), disp(sprintf('\n\n Empty feasible region\n')); return
else, if ~isempty(incomv), ncomv = 1; end; end
A = A(1:m+1,1:nc); A =[A(1:m+1,1:lenA) A(1:m+1,nc)];
[maux, A, z] = compsim(A, maux, 1, 2);
end
if (z == inf | z == -inf), return; end
[m, n] = size(A); x = zeros(n,1); x(maux) = A(1:m-1,n);
x = x(1:n1); z = -A(m,n); t = find(A(m,1:n-1) == 0);
if length(t) > m-1, disp('There are infinite solutions'); end
if ncomv == 1, disp('Redundant constraint(s).'); end
xopt = x; fopt = z;
end

function [maux, A, z]= compsim(A, maux, k, ph)
% Main loop of the simplex primal algorithm.
% Bland's rule to prevente cycling is used.
[m, n] = size(A); [mi, col] = Bland(A(m,1:n-1));
while ~isempty(mi) & mi < 0 & abs(mi) > eps
t = A(1:m-k,col);
if all(t <= 0)
z = -inf; disp(sprintf('\n Unbounded optimal solution with z=%s\n',z));
return
end
[row, small] = minrtest(A(1:m-k,n),A(1:m-k,col));
if ~isempty(row)
if abs(small) <= 100*eps & k == 1, [s,col] = Bland(A(m,1:n-1)); end
A(row,:)= A(row,:)/A(row,col); maux(row) = col;
for i = 1:m, if i ~= row, A(i,:)= A(i,:)-A(i,col)*A(row,:); end; end
[mi, col] = Bland(A(m,1:n-1));
end
end
z = A(m,n);
end

function [m, j] = Bland(D)
% Apply the Bland's rule to the array D.
% m: first negative number in D, j: index of the entry m.
ind = find(D < 0);
if ~isempty(ind), j = ind(1); m = D(j);
else, m = []; j = []; end
end

function [row, mi] = minrtest(a, b)
% Minimum ratio test on vector a and vector b.
% row: index of the pivot row, mi: value of the minimum ratio.
m = length(a); c = 1:m; a = a(:); b = b(:); l = c(b > 0);
[mi, row] = min(a(l)./b(l)); row = l(row);
end
