crit = 1e-6; x1 = -5; h = 0.1; fobj = @(x) 2*x^2*sin(1.5*x) - 4*x^2 + 3*x-1;
[x, f, n] = gsopt(fobj, x1, h, crit)

x1 = 1; [x, f, n] = gsopt(fobj, x1, h, crit)

x1 = 6; [x, f, n] = gsopt(fobj, x1, h, crit)
