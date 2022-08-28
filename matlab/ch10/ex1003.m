x1 = 0.01; h = 0.2; crit = 1e-6; fun = @(x) exp(x) -3*x + 0.02/x - 0.00004/x;
[xopt, fopt, nf] = brentopt(fun, x1, h, crit)
