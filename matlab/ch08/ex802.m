% htcond.m
D = 0.05; h = 98.6; k = 1.18; alpa = 4.97e-7; Ta = 25; Ti = 340; t = linspace(0,15*60,2000);
f = @(x) x*tan(x) - h*D/k/2; x0 = 0.1; gam = fzero(f,x0);
T = Ta + (Ti-Ta)*4*sin(gam).*exp(-4*gam^2*alpa*t/D^2)/(2*gam+sin(2*gam));
fprintf('gamma = %g\n', gam);
plot(t/60,T), grid,xlabel('t(min)'), ylabel('T(deg.C)')
