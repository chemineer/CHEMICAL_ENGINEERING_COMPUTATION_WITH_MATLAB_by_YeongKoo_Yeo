fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2; x0 = [-1.2 1]; crit = 1e-6; kmax = 1e3;
[xopt,fopt, iter] = newtonopt(fun, x0, crit, kmax)
