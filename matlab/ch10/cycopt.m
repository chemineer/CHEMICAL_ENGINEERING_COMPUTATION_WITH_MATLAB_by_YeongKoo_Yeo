function [xopt,fopt,iter] = cycopt(fcyc,x0,crit)
% cycopt.m: minimization by the cyclic coordinate search
% Inputs:
% fcyc: objective functions
% x0: starting point
% crit: stopping criterion
% Outputs:
% xopt: optimal point
% fopt: objective function value at x=xopt
% iter: number of iterations
% Example:
% x0 = [-3 -1 0 1]; crit = 1e-4;
% fc = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
% [xopt,fopt,iter] = cycopt(fc,x0,crit)
stsize = 0.1; x = x0; xk = x; nv = length(x); fk = fcyc(xk); fold = fk; iter = 0;
while (1)
iter = iter + 1;
for j = 1:nv, h(j) = 0; end
for k = 1:nv
for j = 1:nv, d(j) = 0; end
d(k) = 1; aut = 0;
[stfit,fvfit] = quadfit(fcyc,aut,fk,stsize,xk,d,crit);
xk(k) = xk(k) + stfit; h(k) = stfit; fk = fvfit;
end
aut = 0; [stfit,fvfit] = quadfit(fcyc,aut,fk,stsize,xk,h,crit);
for j = 1:nv, xk(j) = xk(j) + stfit * h(j); end
if (abs(fk-fold) <= crit), break; end
fold = fk;
end
xopt = xk; fopt = fcyc(xk);
end

function [x2, f2] = quadfit(fcyc,x1,f1,stsize,x0,d,crit)
fr = 0.05; % 0 < sfrac < 0.5
[x1,x2,x3,f1,f2,f3] = approx3pt(fcyc,x1,f1,stsize,x0,d); % 3-point pattern
tau = (sqrt(5) - 1)/2; redi = 2;
if (x3 < x1), xtemp = x1; ftemp = f1; x1 = x3; f1 = f3; x3 = xtemp; f3 = ftemp; end
iflag = 0; indc = 0; jndc = 0;
while (1)
xdold = abs(x3 - x1); favg = (f1 + f2 + f3)/3.;
if iflag == 0
A = (x1-x2)*(x1-x3); B = (x2-x1)*(x2-x3); C = (x3-x1)*(x3-x2);
x4 = (f1*(x2+x3)/A + f2*(x1+x3)/B + f3*(x1+x2)/C)/(f1/A+f2/B+f3/C)/2;
else
if x2 <= (x1+x3)/2, x4 = x2 + (1-tau)*(x3-x2);
else, x4 = x3 - (1-tau)*(x2-x1);
end
iflag = 0;
end
delt = fr*min(abs(x2-x1),abs(x3-x2));
if abs(x4-x1) < delt, x4 = x1+delt;
elseif abs(x4 - x3) < delt, x4 = x3 - delt;
elseif abs(x4-x2) < delt
if x2 > (x1 + x3)/2, x4 = x2-delt; else, x4 = x2+delt; end
end
f4 = fcyc(x0 + x4*d);
if (x4 > x2)
if (f4 >= f2), x3 = x4; f3 = f4; else, x1 = x2; f1 = f2; x2 = x4; f2 = f4; end
else
if (f4 >= f2), x1 = x4; f1 = f4; else, x3 = x2; f3 = f2; x2 = x4; f2 = f4; end
end
xdnew = abs(x3-x1); fvnew = (f1 + f2 + f3)/3.;
if abs(x3 - x1) <= crit, break; end
if abs(fvnew - favg) <= crit
jndc = jndc + 1; if jndc == 2, break; end
else, jndc = 0;
end

if xdnew/xdold > tau
indc = indc + 1; if indc == redi, indc = 0; iflag = 1; end
else, indc = 0; iflag = 0;
end
end
end
function [x1,x2,x3,f1,f2,f3] = approx3pt(fcyc,x1,f1,stsize,x0,d)
% Three-point pattern
tau = (sqrt(5) - 1)/2; stlen = stsize; x2 = x1 + stlen; x = x0 + x2*d; f2 = fcyc(x);
if (f2 > f1)
temp = x1; x1 = x2; x2 = temp; temp = f1; f1 = f2; f2 = temp; stlen = -stlen;
end
while (1)
stlen = stlen/tau; x3 = x2 + stlen; x = x0 + x3*d; f3 = fcyc(x);
if (f3 > f2), break; else, f1 = f2; x1 = x2; f2 = f3; x2 = x3; end
end
end
