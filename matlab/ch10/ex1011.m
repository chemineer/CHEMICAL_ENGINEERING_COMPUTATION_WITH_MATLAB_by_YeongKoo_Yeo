A = [1 1 1 1 0; 1 2 3 0 1]; b = [7; 12]; c = [-2 -1 -4]; tol = 1e-6;
[xopt,fopt, basic] = barnslp(A,b,c,tol)
