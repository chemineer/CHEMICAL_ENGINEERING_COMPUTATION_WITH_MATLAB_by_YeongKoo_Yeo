xl = -10*[1 1]; xu = 10*[1 1]; nb = 8; ps = 50; ng = 60; mp = 0.05;
[xopt,fopt,iter] = gaopt(@fcy2,xl,xu,nb,ps,ng,mp)
