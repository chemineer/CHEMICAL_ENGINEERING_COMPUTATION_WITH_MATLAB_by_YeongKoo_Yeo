function Pv = vpwagner(T,Tc,Pc,C) 
% Estimation of vapor pressure using the Wagner equation 
% input: 
% T,Tc: temperature (K) and critical temperature (K) 
% Pc: critical pressure (MPa) 
% C: parameter vector of the Wagner equation 
% output 
% Pv: estimated vapor pressure (MPa) 
Tr = T./Tc; a = C(1); b = C(2); c = C(3); d = C(4); 
Pv = Pc.*exp((a*(1-Tr) + b*(1-Tr).^1.5 + c*(1-Tr).^2.5 + d*(1-Tr).^5)/Tr); 
fprintf('Vapor pressure = %g MPa \n', Pv); 
end 