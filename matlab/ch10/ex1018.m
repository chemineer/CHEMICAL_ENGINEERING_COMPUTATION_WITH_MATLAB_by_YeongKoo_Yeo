f = @(x) 100*(x(1)^2 - x(2))^2 + (1 - x(1))^ 2; x0 = [-1 1]; crit = 1e-6;
[xopt,fopt,iter] = hjopt(f,x0,crit)

