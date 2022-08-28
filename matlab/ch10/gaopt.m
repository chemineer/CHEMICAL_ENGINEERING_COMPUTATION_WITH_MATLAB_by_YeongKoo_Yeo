function [xopt,fopt,iter] = gaopt(fun,xl,xu,nb,ps,ng,mp)
% gaopt.m: maximization by the genetic algorithm
% Inputs:
% fun: objective functions
% xl,xu: lower and upper bounds on x
% nb: number of binary digits
% ps: population size
% ng: number of generations
% mp: mutation probability
% Outputs:
% xopt: optimal point
% fopt: objective function value at x=xopt
% iter: number of function evaluations
% Example:
% xl=-6*[1 1 1 1]; xu=6*[1 1 1 1]; nb=8; ps=40; ng=50; mp=0.02;
% f = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
% [xopt,fopt,iter] = gaopt(f,xl,xu,nb,ps,ng,mp)
n = length(xu);
Pn = round(rand(ps,n*nb)); % initial population
for kg = 1:ng
ft = fitness(fun,ps,nb,n,xl,xu,Pn); [fmax,kmax] = max(ft);
if (kg == 1)
fbest = fmax;
for k = 1:n
ab = num2str(Pn(kmax,nb*(k-1)+1:nb*k));
adec = bin2dec(ab); % convert binary string to decimal integer
xbest(k) = xl(k) + (xu(k)-xl(k))/(2^nb - 1)*adec;
end
else
if (fmax > fbest)
fbest = fmax;
for k = 1:n
ab = num2str(Pn(kmax,nb*(k-1)+1:nb*k)); adec = bin2dec(ab);
xbest(k) = xl(k) + (xu(k)-xl(k))/(2^nb - 1)*adec;
end
end
end
[frk,ft] = frank(ps, ft); % scaling
for k = 1: ps % shuffling: roulette wheel
ik = frk(k); ptemp = ps; ktemp = k; am(ik) = 2*(ptemp + 1 - ktemp)/ (ptemp*(ptemp + 1));
end
for k = 2: ps, ptemp = am(k - 1); am(k) = ptemp + am(k); end
%
for k = 1:ps % selection
kstr1 = roul(ps,am); kstr2 = roul(ps,am); kchd = cros(Pn,nb,n,kstr1, kstr2);
kchd = mutat(kchd,mp,nb,n);
for j = 1: nb*n, Pn(k, j) = kchd(j); end
end
end
xopt = xbest; fopt = fbest; iter = ng*ps;
end
function ft = fitness(fun,ps,nb,n,xl,xu,Pn) % function evaluation
for k = 1:ps
for j = 1:n
a = num2str(Pn(k,nb*(j-1) + 1:nb*j)); adec = bin2dec(a);
x(j) = xl(j) + (xu(j)-xl(j))/(2^nb - 1)*adec;
end
f = fun(x); ft(k) = f;
end
end

function [frk,ft] = frank(ps, ft)
for k = 1:ps, frk(k) = k; end
for k = 1:ps-1
temp = ft(k); ktemp = k;
for j = k+1:ps
if (ft(j) > temp), temp = ft(j); ktemp = j; end
end
jtemp = frk(ktemp); frk(ktemp) = frk(k); frk(k) = jtemp;
ftemp = ft(ktemp); ft(ktemp) = ft(k); ft(k) = ftemp;
end
end
function kstr = roul(ps,am) % shuffling operation
temp = rand;
for k = 1:ps
while (1)
if (temp > am(k)), break; end
kstr = k; return;
end
end
end
function kchd = cros(Pn,nb,n,kstr1, kstr2) % crossover operation
for j = 1:n
jp = nb*(j - 1); nc = floor(rand*nb + 0.5);
if (nc == 0)
for k = 1:nb, kchd(k+jp) = Pn(kstr2, k+jp); end
elseif (nc == nb)
for k = 1:nb, kchd(k+jp) = Pn(kstr1, k+jp); end
else
for k = 1: nc, kchd(k+jp) = Pn(kstr1, k+jp); end
for k = nc+1:nb, kchd(k+jp) = Pn(kstr2, k+jp); end
end
end
end
function [kchd] = mutat(kchd,mp,nb,n) % mutation operation
for k = 1:nb*n
if (rand <= mp)
if (kchd(k) == 1), kchd(k) = 0;
else, kchd(k) = 1;
end
end
end
end
