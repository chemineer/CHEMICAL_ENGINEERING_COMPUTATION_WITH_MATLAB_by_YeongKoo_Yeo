% membrx.m
Fa0 = 10; k = 0.7; kc = 0.2; Kc = 0.05; R = 8.314; P = 830.6; T = 500;
Ct0 = P/(R*T); Vf = 500; Vc = 400; x0 = [Fa0 0 0];
[V x] = ode45(@mrf,[0 Vf],x0,[],k,kc,Kc,Ct0);
Fa = x(:,1); Fb = x(:,2); Fc = x(:,3);
fprintf('Final values: Fa = %g, Fb = %g, Fc = %g\n', Fa(end), Fb(end), Fc(end));
plot(V,Fa,V,Fb,':',V,Fc,'--'), xlabel('V(dm^3)'), axis tight
ylabel('Molar flow rate (mol/min)'), legend('F_A','F_B','F_C')
for i = 1:length(V)
if V(i) >= Vc, iV = i; break; end
end
Fc = Fa(iV) + (Fa(iV+1)-Fa(iV))*(Vc - V(iV))/(V(iV+1)-V(iV)); Xa = (Fa0 - Fc)/Fa0;
fprintf('At V = %g, Fa = %g and the conversion of A = %g\n', Vc, Fc, Xa);
function dxdv = mrf(v,x,k,kc,Kc,Ct0)
% x(1) = Fa, x(2) = Fb, x(3) = Fc
Ft = sum(x); rA = -k*Ct0*(x(1)/Ft - Ct0*x(2)*x(3)/(Kc*Ft^2));
dxdv = [rA; -rA - kc*Ct0*x(2)/Ft; -rA];
end
