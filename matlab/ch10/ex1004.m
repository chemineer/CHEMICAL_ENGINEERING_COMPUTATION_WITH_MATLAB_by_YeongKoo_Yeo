C = 8; a = -3; b = 8; crit = 1e-6; nfmax = 2000; fun = @(x) -sin(1.2*x)- sin(3.5*x);
[xopt,fopt,nf] = shubertopt(fun, a, b, C, crit, nfmax)
