fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
delfun = @(x) [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1); 200*(x(2)-x(1)^2)];
x0 = [-1.2 1]; crit = 1e-6; alpha0 = 5; kmax = 1e4;
[xopt,fopt, iter] = sdopt(fun, delfun, x0, alpha0, crit, kmax)
