fun = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2+90*(x(4)-x(3)^2)^2+(1-x(3))^2+10.1*((x(2)-1)^2+(x(4)-1)^2)+19.8*(x(2)-1)*(x(4)-1);
delfun = @(x) [400*x(1)*(x(1)^2-x(2))+2*(x(1)-1); 200*(x(2)-x(1)^2)+20.2*(x(2)-1)+19.8*(x(4)-1); 360*x(3)*(x(3)^2-x(4))+2*(x(3)-1); 180*(x(4)-x(3)^2)+20.2*(x(4)-1)+19.8*(x(2)-1)];
x0 = [-3 -1 -3 -1]; crit = 1e-6; alpha0 = 2; kmax = 1e3;
[xopt,fopt, iter] = dfpopt(fun, delfun, x0, alpha0, crit, kmax)