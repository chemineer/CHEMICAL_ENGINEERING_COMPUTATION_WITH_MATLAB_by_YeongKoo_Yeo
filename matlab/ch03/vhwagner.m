function dHv = vhwagner(T,Tc,Pc,vL,vV,C) 
% Estimation of enthalpy of vaporization using the Wagner equation 
% input: 
% T,Tc: temperature and critical temperature (K) 
% Pc: critical pressure (MPa) 
% C: parameter vector of the Wagner equation 
% output 
% dHv: estimated enthalpy of vaporization (J/mol) 
Tr = T./Tc; a = C(1); b = C(2); c = C(3); d = C(4); 
Pv = Pc.*exp((a*(1-Tr) + b*(1-Tr).^1.5 + c*(1-Tr).^2.5 + d*(1-Tr).^5)/Tr); 
wd = log(Pv./Pc) + a + 1.5*b*(1-Tr).^0.5 + 2.5*c*(1-Tr).^1.5 + 5*d*(1-Tr).^4; 
dHv = -Pv*(vV - vL).*wd*1e6; % J/mol 
end 