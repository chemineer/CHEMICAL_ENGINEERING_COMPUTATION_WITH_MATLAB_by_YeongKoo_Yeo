fun = @(x) 0.02*x(1)^2+1.2*x(2)^2-80; delfun = @(x) [0.04*x(1); 2.4*x(2)];
ne = 0; m = 6; crit = 1e-4; x0 = [4 5]; A = [-1 0;-10 1;-1 0;0 -1;1 0;0 1]; b = [-3 -12 40 40 40 40]';
[xopt,fopt,iter] = rosencopt(fun,delfun,x0,A,b,ne,m,crit)
