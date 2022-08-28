x0 = [2 3]'; lam0 = 1; mu0 = 3; crit = 1e-6;
[xopt,fopt,iter] = sqpopt(@fun,@dfun,x0,lam0,mu0,crit)
