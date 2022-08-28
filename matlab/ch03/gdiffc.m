function Df = gdiffc(T,M,phi,mu,v) 
% Estimation of diffusion coefficients in liquid phase 
% input 
% T: temperature(C) (scalar or row vector) 
% M: molecular weight
% phi association factor of solvent B 
% mu: viscosity of solvent B (cP) 
% v: molal volume of solute A at its normal boiling point (cm^3/gmol) 
% output 
% Df: estimated diffusion coefficient in liquid phase (cm^2/s) 
Df = 7.4e-8*T.* sqrt(phi.*M)./(mu.*v.^0.6); 
end 