% grgoptex.m
x0 = [0 0 0 0 7 11 6]; xl = [-100 -100 -100 -100 0 0 0]; xu = 100*ones(1,7); kmax = 1e3; crit = 1e-4;
[xopt,fopt,iter] = grgopt(@grgfun,@delgrgf,@delgrgg,x0,xl,xu,kmax,crit)

function df = delgrgf(x)
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); df = zeros(1,length(x));
df(1) = 2.3*x1-4; df(2) = 2*x2-6; df(3) = 4.6*x3-20; df(4) = -2.4*x4+6;
end

function dg = delgrgg(x)
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
dg = [2*x1+1 2*x2-1 3*x3+1 2*x4-1 1 0 0;
2*x1-1 4*x2 2*x3 4*x4-1 0 1 0;
4*x1+2 2*x2-1 2*x3 -1 0 0 1];
end
