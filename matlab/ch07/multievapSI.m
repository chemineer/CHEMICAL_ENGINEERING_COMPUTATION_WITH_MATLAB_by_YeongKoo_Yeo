function res = multievapSI(evdat)
% Multiple-effect evaporator calculation (SI unit system)
% Data
crit = 1e-3; Tc = 647.096; % critical temperature of H2O (K)
Pc = 22064000; % critical pressure of H2O (Pa)
T0 = 273.15; xf = evdat.xf; xpn = evdat.xp; mf = evdat.mf;
Tf = evdat.Tf + T0; % C -> K
Ps = evdat.Ps; % steam pressure (Pa)
Pn = evdat.Pn; % pressure of the final effect (Pa)
U = evdat.U; % overall heat transfer coefficient
n = length(U);
% Saturation temperature
a = [-7.85951783 1.84408259 -11.7866497 22.6807411 -15.9618719 1.80122502];
gs = @(x) a(1)*x + a(2)*x^1.5 + a(3)*x^3 + a(4)*x^3.5 + a(5)*x^4 + a(6)*x^7.5 - (1-x)*log(Ps/Pc);
gn = @(x) a(1)*x + a(2)*x^1.5 + a(3)*x^3 + a(4)*x^3.5 + a(5)*x^4 + a(6)*x^7.5 - (1-x)*log(Pn/Pc);
xs = fzero(gs,0.5); xn = fzero(gn,0.5);
Ts = Tc*(1 - xs); % steam temperature (K)
Tn = Tc*(1 - xn); % temperature of the final effect (K)
% Enthalpy
sts = satsteam(Ts-T0); % saturated steam
stv = satsteam(Tn-T0); % final effect
stf = satsteam(Tf-T0); % feed
Hv0 = sts.hV; hp0 = sts.hL; % J/kg
Hvn = stv.hV; hpn = stv.hL; hf = stf.hL;
% Mass balances and temperature drop
mpn = xf*mf/xpn; % kg/hr
dTtotal = Ts - Tn; % K
sumU = sum(1./U);
dT = (1./U)*dTtotal/sumU; % K
critA = 10; oldA = 10*ones(1,n); iter = 0;
while critA >= crit
T(1) = Ts - dT(1);
for j = 2:n, T(j) = T(j-1) - dT(j); end
for j = 1:n
sprop(j) = satsteam(T(j)-T0); Hv(j) = sprop(j).hV; hp(j) = sprop(j).hL;
end
% balance equations
evM = [Hv0-hp0 hp(1)-Hv(1) 0; 0 Hv(1)+hp(2)-2*hp(1) hp(2)-Hv(2);...
0 Hv(3)-hp(2) Hv(2)+Hv(3)-2*hp(2)];
evb = [hp(1)-hf; hp(2)-hp(1); Hv(3)-hp(2)+(hp(3)-Hv(3))*xf/xpn]*mf;
mv = evM\evb; q(1) = mv(1)*(Hv0-hp0);
for j = 2:n, q(j) = mv(j)*(Hv(j)-hp(j)); end
Area = q./(U.*dT); avgA = sum(Area)/n; critA = sum(abs(Area - oldA)); dT = dT.*Area/avgA;
if abs(sum(dT) - dTtotal) >= crit, dT = dT*dTtotal/sum(dT); end
iter = iter + 1; oldA = Area;
end
% Results
res.T = T-T0; % vector of temperature of each effect (C)
res.A = Area; % vector of area of each effect (m^2)
res.mv = mv; % vector of vapor flow rate from each effect (kg/hr)
res.iter = iter; % number of iterations
end
