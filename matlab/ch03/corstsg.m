function sg = corstsg(Pc,Tc,Tb,T) 
% Estimation of surface tension using the corresponding states correlation 
% input 
% Pc: critical pressure (bar) 
% Tc: critical temperature (K) 
% Tb: normal boiling point (K) 
% T: temperature (K) 
% output 
% sg: surface tension (dyne/cm) 
Tr = T./Tc; Tbr = Tb./Tc; Q = 0.1196*(1 + Tbr.*log(Pc/1.01325)./(1-Tbr)) - 0.279; 
sg = Pc.^(2/3).*Tc.^(1/3).*Q.*(1-Tr).^(11/9); 
fprintf('Surface tension = %g dyne/cm \n', sg); 
end 