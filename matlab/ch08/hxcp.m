function Cp = hxcp(cpref,Tref,T)
% Calculate heat capacity at T using two ref. temperatures and heat capacities
% input:
%  cpref: heat capacity vector (J/kg/K) at ref. temperature vector Tref
%  Tref: reference temperature vector (K)
%  T: temperature (K) at which heat capacity is to be determined
% output:
%  Cp: heat capacity (J/kg/K) (Cp = A+B*T)
Cp = (cpref(2)*Tref(1) - cpref(1)*Tref(2) + T*(cpref(1)-cpref(2)))/(Tref(1) - Tref(2));
end
