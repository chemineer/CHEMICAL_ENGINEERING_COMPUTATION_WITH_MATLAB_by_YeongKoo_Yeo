function mu = hxvis(muref,Tref,T,fstate)
% Calculate viscosity at T using two reference temperatures and viscosities
% input:
%  muref: viscosity vector at reference temperature vector Tref
%  Tref: reference temperature vector
%  T: temperature (K) at which viscosity is to be determined
%  fstate: state of fluid (1: liquid(mu = A*exp(B/T), 2: vapor(mu = A+BT))
% output:
%  mu: viscosity (Ns/m^2)
if fstate == 1 % liquid
mu = muref(1)*exp(Tref(2)*(Tref(1) - T)/(T*Tref(1) - Tref(2))*log(muref(2)/ muref(1)));
else % vapor
mu = (muref(2)*Tref(1) - muref(1)*Tref(2) + T*(muref(1)-muref(2)))/(Tref(1) - Tref(2));
end
end
