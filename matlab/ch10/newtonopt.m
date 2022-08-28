function [xopt,fopt, iter] = newtonopt(fun, x0, crit, kmax)
% newtonopt.m: minimization by Newton's method
% inputs:
% fun: objective function
% x0: starting point
% crit: stopping criterion
% kmax: maximum iterations
% outputs:
% xopt: optimal point
% fopt: function value at the optimal point (=f(xopt))
% iter: number of iterations
% Example:
% fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% x0 = [-1.2 1]; crit = 1e-6; kmax = 1e3;
% [xopt,fopt, iter] = newtonopt(fun, x0, crit, kmax)

h = 1e-4; fx = fun(x0); nf = length(fx); nx = length(x0);
xs(1,:) = x0(:)';  % initial row solution vector
for k = 1:kmax
dx = -jacob(fun, xs(k,:), h)\fx(:); % -[df]^(-1)*fx
xs(k+1,:) = xs(k,:) + dx'; fx = fun(xs(k+1,:));
if norm(fx) < crit | norm(dx) < crit, break; end
end
x = xs(k+1,:); xopt = x; fopt = fx; iter = k;
end


