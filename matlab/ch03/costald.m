function rhoL = costald(T,Tc,vc,w,Mw) 
% estimation of liquid density by COSTALD method 
% input: 
% T,Tc: temperature and critical temperature (K) 
% vc: critical volume 
% w: acentric factor 
% Mw: molecular weight (g/mol) 
% output 
% rhoL: estimated density (kg/m^3) 
a = -1.52816; b = 1.43907; c = -0.81446; d = 0.190454;  
e = -0.296123; f = 0.386914; g = -0.0427458; h = -0.0480645; 
Tr = T./Tc; Vr0 = 1+a.*(1-Tr).^(1/3)+b.*(1-Tr).^(2/3)+c.*(1-Tr)+d.*(1-Tr).^(4/3); 
Vrd = (e + f*Tr + g*Tr.^2 + h*Tr.^3)./(Tr - 1.00001); 
spV = vc.*Vr0.*(1 - w.*Vrd); % specific volume(cm^3/mol) 
rhoL = 1000*Mw/spV; % density(kg/m^3) 
fprintf('Density = %g g/cm^3\n', rhoL) 
end 