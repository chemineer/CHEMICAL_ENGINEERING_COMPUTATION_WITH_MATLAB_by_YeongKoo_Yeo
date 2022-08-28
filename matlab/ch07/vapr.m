function dz = vapr(t,z)
% define differential equations: z(1)=mV, z(2)=mL
vapdat;  % retrieve data
mV = z(1); mL = z(2); Vv = V - mL/rhol;
% Determine T by using bisection method
f = @(T) exp(A-B/(T+C)) - 7.6e5*mV*R*T/(Vv*Mw); Ta = 273.15; Tb = Ta + 400;
if f(Ta)*f(Tb) > 0, display('No solution T'), return; end
Tm = (Ta + Tb)/2;
while abs(Ta-Tb) > 1e-2
if f(Ta)*f(Tm) < 0, Tb = Tm; else, Ta = Tm; end; Tm = (Ta + Tb)/2;
end
T = Tm;
% Define differential eqns
P = exp(A - B/(T + C)) / 750.044; % P: bar
vB = (Fi*rhol*Cp*Ti + Q)/(lambda + Cp*(T - 273.15));
vO = Kv*sqrt(P*(P - P0));
dz = [vB - vO;  Fi*rhol - vB];
end
