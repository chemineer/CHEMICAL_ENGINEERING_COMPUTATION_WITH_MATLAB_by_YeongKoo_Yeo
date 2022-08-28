f = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
x0 = [-3 -1 0 1]; crit = 1e-6; [xopt,fopt,iter] = nmopt(f,x0,crit)
