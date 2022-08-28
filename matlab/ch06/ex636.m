% geferm.m
clear all;
mumax=0.33; Cps=93; Ks=1.7; kd=0.01; m=0.03; Ysc=1/0.08; Ypc=5.6; % data
C0 = [1 250 0]; tspan = [0 12]; % initial concentration and time span
[t C] = ode45(@gf,tspan,C0,[],mumax,Cps,Ks,kd,m,Ysc,Ypc); Cc = C(:,1); Cs = C(:,2); Cp = C(:,3);
rd = kd*Cc; rsm = m*Cc; rg = mumax*(1-Cp/Cps).^0.52 .* Cc.*Cs./(Ks+Cs);
subplot(2,2,1), plot(t,Cc), xlabel('t(hr)'), ylabel('C_c(g/dm^3)')
subplot(2,2,2), plot(t,Cs,t,Cp,'--'), xlabel('t(hr)'), ylabel('C(g/dm^3)'),
legend('C_s','C_p','Location','best')
subplot(2,2,3), plot(t,rg,t,rsm,'.-',t,rd,'--'), xlabel('t(hr)'), ylabel('rates(g/dm^3/hr)')
legend('r_g','r_sm','r_d','Location','best')
fprintf('Final concentrations: Ccf = %g, Csf = %g, Cpf = %g\n',Cc(end),Cs(end),Cp(end));
fprintf('Final reaction rates: rgf = %g, rsmf = %g, rdf = %g\n',rg(end),rsm(end),rd(end));
function dC = gf (t,C,mumax,Cps,Ks,kd,m,Ysc,Ypc)
% C(1)=Cc, C(2)=Cs, C(3)=Cp
rd = kd*C(1); rsm = m*C(1); rg = mumax*(1 - C(3)/Cps)^0.52 * C(1)*C(2)/(Ks + C(2));
dC = [rg - rd; -rg*Ysc - rsm; Ypc*rg];
end
