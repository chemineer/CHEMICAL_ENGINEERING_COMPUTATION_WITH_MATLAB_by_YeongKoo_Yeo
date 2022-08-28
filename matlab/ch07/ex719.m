% compdyndist.m
dpar.alpha = 1.5; dpar.n = 30; dpar.nf = 15; dpar.F = 1; dpar.zf = 0.5; dpar.q = 1;
dpar.R = 2.7;
dpar.Vs = 3.2; dpar.D = 0.5; dpar.md = 5; dpar.mb = 5; dpar.mt = 0.5; dels.delR = 0.01*dpar.R;
dels.delRt = 10; dels.delV = 0; dels.delVt = 0; dels.delz = 0; dels.delzt = 0;
dels.delF = 0;
dels.delFt = 0; t0 = 0; tf = 400; x0 = 0.5*ones(1,dpar.n); nv = 1:dpar.n;
x0 = fsolve(@ssdist,x0,[],dpar); % steady-state to be used as initial conditions
for i = 1:length(x0), if x0(i) <= 0, x0(i) = -x0(i); end; end
[t x] = ode45(@dyndist,[t0 tf],x0,[],dpar,dels);
subplot(1,2,1), plot(t,x(:,1),t,x(:,end),'--'), xlabel('t(min)'), ylabel('x'),legend('x_D','x_B')
subplot(1,2,2), plot(nv,x(end,:)), xlabel('n'), ylabel('x_i'), axis tight
