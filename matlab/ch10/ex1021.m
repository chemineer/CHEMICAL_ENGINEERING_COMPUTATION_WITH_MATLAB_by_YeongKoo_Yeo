T = 100; xl = -10*[1 1]; xu = 10*[1 1]; rp = 0.8; rs = 0.9; crit = 1e-8;
[xopt,fopt,iter] = saopt(@fcy,T,xl,xu,rp,rs,crit)
