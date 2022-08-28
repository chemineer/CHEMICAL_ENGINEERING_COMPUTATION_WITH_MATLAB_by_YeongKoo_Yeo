% microrx.m
Fa0 = 2.26e-5; R = 8.314; P0 = 1641; T0 = 698; E = 24000;
Ct0 = P0/(R*T0); Vf = 1e-5; x0 = [Fa0 0 0];
k = 0.29*exp(E/1.987*(1/500-1/T0));
dxdv = @(t,x) [-k*Ct0^2*(x(1)/ sum(x))^2; k*Ct0^2*(x(1)/ sum(x))^2; (k/2)*Ct0^2*(x(1)/ sum(x))^2];
[V x] = ode23s(dxdv,[0 Vf],x0); Fa = x(:,1); Fb = x(:,2); Fc = x(:,3);
fprintf('Final values: Fa = %g, Fb = %g, Fc = %g\n', Fa(end), Fb(end), Fc(end));
plot(V,Fa,V,Fb,':',V,Fc,'--'), xlabel('V(dm^3)'), axis tight
ylabel('Molar flow rate (mol/min)'), legend('F_A','F_B','F_C')
