function dT = pflowht(x,T,pf)
% heat transfer in fluid flowing through a pipe
% Retrieve data
r = pf.r; v = pf.v; n = pf.n; h = pf.h; alpa = pf.alpa; Tb = pf.Tb; dT = zeros(n,1);
% Difference model equations
for k = 1:n
if k == 1, s = 2*(T(k+1)-T(k))/h^2; d = 0;
elseif k == n, s = (Tb- 2*T(k)+T(k-1))/h^2; d = (Tb-T(k-1))/(2*h*r(k+1));
else, s = (T(k+1)-2*T(k)+T(k-1))/h^2; d = (T(k+1)-T(k-1))/(2*h*r(k)); end
dT(k) = (alpa/v(k))*(s + d);
end
end
