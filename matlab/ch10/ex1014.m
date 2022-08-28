% optzout
delf = @(x) [-1.2; -3]; delg = @(x) [2*x(1); 12*x(2)];
x0 = [1 0]; xl = [0 0]; xu = [10 10]; nc = 0; ncs = 1; crit = 1e-4; kmax = 1e3;
[xopt,fopt,iter] = zoutopt(@funz,delf,delg,x0,xl,xu,nc,ncs,crit,kmax)

function fg = funz(x)
fg(1) = -(1.2*x(1) + 3*x(2));
fg(2) = x(1)^2 + 6*x(2)^2 - 1;
end
