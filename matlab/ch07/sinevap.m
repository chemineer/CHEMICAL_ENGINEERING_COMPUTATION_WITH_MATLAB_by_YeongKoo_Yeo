function res = sinevap(evdat)
% Calculation of single-effect evaporator
% Data
Tc = 647.096; % critical temperature of H2O (K)
Pc = 22064000; % critical pressure of H2O (Pa)
T0 = 273.15; xf = evdat.xf; xp = evdat.xp; mf = evdat.mf;
Tf = (evdat.Tf-32)/1.8 + 273.15; % F->K
Ps = evdat.Ps*6894.757; Pv = evdat.Pv*6894.757; % psia->Pa
U = evdat.U; cf = 1/2326; % J/kg->Btu/lb
% Saturation temperature
a = [-7.85951783 1.84408259 -11.7866497 22.6807411 -15.9618719 1.80122502];
gs = @(x) a(1)*x + a(2)*x^1.5 + a(3)*x^3 + a(4)*x^3.5 + a(5)*x^4 +a(6)*x^7.5 - (1-x)*log(Ps/Pc);
gv = @(x) a(1)*x + a(2)*x^1.5 + a(3)*x^3 + a(4)*x^3.5 + a(5)*x^4 +a(6)*x^7.5 - (1-x)*log(Pv/Pc);
xs = fzero(gs,0.5); xv = fzero(gv,0.5); Ts = Tc*(1 - xs); Tv = Tc*(1 - xv);
% Enthalpy
sts = satsteam(Ts-T0); stv = satsteam(Tv-T0); stf = satsteam(Tf-T0);
Hs = sts.hV*cf; hs = sts.hL*cf; % steam
Hv = stv.hV*cf; hv = stv.hL*cf; % solution
hf = stf.hL*cf;
% Material balance
mp = xf*mf/xp; % Btu/hr
mev = mf - mp; q = mev*Hv - mf*hf + mp*hv; ms = q/(Hs - hs);
Ts = (Ts - 273.15)*1.8+32; % K->F
Tv = (Tv - 273.15)*1.8+32; dT = Ts - Tv; A = q/U/dT;
% Assign output structure (heating surface(A) and steam amount(ms))
res.ms = ms; res.A = A;
end
