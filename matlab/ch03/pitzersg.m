function sg = pitzersg(w,Pc,Tc,T) 
% Estimation of surface tension using Pitzer's relation 
% input 
% Pc: critical pressure (bar) 
% Tc: critical temperature (K) 
% w: acentric factor 
% T: temperature (K) 
% output 
% sg: surface tension (dyne/cm)  

Tr = T./Tc; 
sg = (Pc.^(2/3)).*(Tc.^(1/3)).*((1.86+1.18*w)/19.05).*((3.75+...    
0.91*w)./(0.291-0.08*w)).^(2/3).*(1-Tr).^(11/9); 
fprintf('Surface tension = %g dyne/cm \n', sg); 
end 