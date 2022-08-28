% stcstr.m 
Ca0 = 0.4; k0 = 460; E = 1380; tau = 0.18; Tc = 298.15; kappa = 78; Cp = 32;  
x0 = [0.1 300]; x = fsolve(@cstrst, x0, [], Ca0, k0, E, tau, Tc, kappa, Cp); 
x(2) = x(2) - 273.15  % K -> deg.C  

function y = cstrst(x,Ca0, k0, E, tau, Tc, kappa, Cp) 
% x(1) = Ca, x(2) = T 
ra = -k0*exp(-E/x(2))*x(1); dHr = -151080 + 2*(x(2) - 298.15); 
y = [(Ca0 - x(1))/tau - (-ra); (-dHr/Cp)*(-ra/Ca0) - (1+kappa)*(x(2) - Tc)/ tau ]; 
end 