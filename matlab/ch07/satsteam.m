function stpr = satsteam(Ts)
% Properties of the saturated steam at a given temperature Ts(C)
% stpr: a structure containing results of calculations
% Critical properties
Tc = 647.096;  %  critical temperature(K)
Pc = 22064000;  % critical pressure(Pa)
rhoc = 322;   %  critical density(kg/m^3)
% Definition of parameters and constants
Ts = Ts + 273.15; alpha0 = 1000; % alpha_0(J/kg)
phi0 = 1000/647.096;
a = [-7.85951783 1.84408259 -11.7866497 22.6807411 -15.9618719 1.80122502];
b = [1.99274064 1.09965342 -0.510839303 -1.75493479 -45.5170352 -674694.45];
c = [-2.03150240 -2.68302940 -5.38626492 -17.2991605 -44.7586581 -63.9201063];
d = [-5.65134998e-8 2690.66631 127.287297 -135.003439 0.981825814];
alphad = -1135.905627715; phid = 2319.5246; theta = Ts/Tc; tau = 1 - theta;
% saturated steam pressure
tw = (Tc/Ts)*(a(1)*tau + a(2)*tau^1.5 + a(3)*tau^3+a(4)*tau^3.5 + a(5)*tau^4 + a(6)*tau^7.5);
Ps = Pc*exp(tw);
% density of saturated liquid
rhoL =rhoc*(1 + b(1)*tau^(1/3) + b(2)*tau^(2/3) + b(3)*tau^(5/3) +...
b(4)*tau^(16/3) + b(5)*tau^(43/3) + b(6)*tau^(110/3));
% density of saturated steam
vw = c(1)*tau^(1/3) + c(2)*tau^(2/3) + c(3)*tau^(4/3) + c(4)*tau^3 + c(5)*tau^(37/6) + c(6)*tau^(71/6);
rhoV = rhoc*exp(vw);
% specific volume
vL = 1/rhoL;  %  saturated liquid
vV = 1/rhoV;  %  saturated steam
% alpha
alpha = alpha0*(alphad + d(1)*theta^(-19) + d(2)*theta + d(3)*theta^4.5 + d(4)*theta^5 + d(5)*theta^54.5);
% phi
phi = phi0*(phid + (19/20)*d(1)*theta^(-20) + d(2)*log(theta) +...
(9/7)*d(3)*theta^3.5 + (5/4)*d(4)*theta^4 + (109/117)*d(5)*theta^53.5);
% dp/dT
tv = 7.5*a(6)*tau^6.5 + 4*a(5)*tau^3 + 3.5*a(4)*tau^2.5 +...
3*a(3)*tau^2 + 1.5*a(2)*tau^0.5 + a(1) + log(Ps/Pc);
dpdT = (-Ps/Ts)*tv;
% enthalpy
hL = alpha + (Ts/rhoL)*dpdT; % saturated liquid
hV = alpha + (Ts/rhoV)*dpdT; % saturated steam
% entropy
sL = phi + (1/rhoL)*dpdT;  % saturated liquid
sV = phi + (1/rhoV)*dpdT; % saturated steam
% result structure
stpr.T = Ts; stpr.P = Ps; stpr.vL = vL; stpr.vV = vV; stpr.hL = hL;
stpr.hV = hV; stpr.sL = sL; stpr.sV = sV;
end
