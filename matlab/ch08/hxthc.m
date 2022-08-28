function xk = hxthc(xkref,Tref,T)
% Calculate heat conductivity at T using two ref. temperatures and
% conductivities
% input:
%  xkref: heat conductivity vector at ref. temperature vector Tref
%  Tref: reference temperature vector (K)
%  T: temperature (K) at which heat conductivity is to be determined
% output:
%  xk: heat conductivity (W/m/K) (xk = A+B*T)
xk = (xkref(2)*Tref(1) - xkref(1)*Tref(2) + T*(xkref(1)-xkref(2)))/(Tref(1) - Tref(2));
end
