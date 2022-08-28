f = @(x) x(1)^2 + 2*x(2)^2 + 2*x(1)*x(2); x0 = [0.5 1]; crit = 1e-6;
[xopt,fopt,istag] = rbopt(f,x0,crit)
