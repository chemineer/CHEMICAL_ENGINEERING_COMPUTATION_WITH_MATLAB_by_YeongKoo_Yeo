A = [1 3 2 4;1 1 1 0;0 0 1 1]; b = [8 2 1]; c = [12 10 6 4]; nl = 1; ne = 2; ng = 0; ibd = [1 2 3 4];
[xopt,fopt,iter] = bnbopt(A,b,c,nl,ng,ne,ibd)
