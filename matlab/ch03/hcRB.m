function cpL = hcRB(T,Tc,w,Cpi) 
% Estimation of liquid heat capacity by Rowlinson/Bondi method 
% input: 
% T,Tc: temperature and critical temperature (K) 
% w: acentric factor 
% Cpi: ideal gas heat capacity (J/mol/K) 
% output 
% cpL: estimated liquid heat capacity (J/mol/K) 
Tr = T./Tc; R = 8.3143; 
cpL = Cpi + 1.45*R + 0.45*R./(1-Tr) + 0.25*w*R*(17.11 +...    
25.2*(1-Tr).^1/3 ./ Tr + 1.742./(1-Tr)); % (J/mol/K) 
fprintf(' Heat capacity = %g J/mol/K \n', cpL) 
end 