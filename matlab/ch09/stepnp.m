function [y,t] = stepnp(num, den, t0, delt, fint, ms);
% Calculate step response of a proper SISO system
% num : numerator of transfer function
% den : denominator of transfer function
% t0 : time at which unit step input is introduced
% delt : time step
% fint : final response time
% y : step response
% ms: step size
% Partial fraction of (transfer function)*(step input)
% (r: residue vector, p: pole vector, k: constant vector)
[r, p, k] = residue(num, conv(den, [1 0]));
% Set calculation time interval
t = t0 : delt : fint;
% Identify pole multiplicity
for j = 1 : size(p)
n = 1;
for i = 1 : size(p)
if p(j) == p(i), if (i ~= j), n = n+1; end; end
end
mult(:, j) = n;
end
% Step response: use inverse Laplace transform
y = zeros(size(t));
j = 1;
while j <= size(p, 1)
for i = 1 : mult(:, j), y = y + r(j+i-1)*((t-t0).^(i-1)).*exp(p(j)*(t-t0))/ factorial(i-1); end
j = j + i;
end
y = ms*y;
end
