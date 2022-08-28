function rho = hxrho(rhoref,Tref,T,fstate)
% Calculate density at T using two ref. temperatures and densities
% input:
%  rhoref: density vector (kg/m^3) at ref. temperature vector Tref
%  Tref: reference temperature vector (K)
%  T: temperature (K) at which density is to be determined
%  fstate: state of fluid (1: liquid(rho = A+B*T, 2: vapor(rho = A/T))
% output:
%  rho: density (kg/m^3)
if fstate == 1 % liquid
rho = (rhoref(2)*Tref(1) - rhoref(1)*Tref(2) + T*(rhoref(1)-rhoref(2)))/ (Tref(1) - Tref(2));
else % vapor
rho = rhoref(1)*Tref(1)/T;
end
end
