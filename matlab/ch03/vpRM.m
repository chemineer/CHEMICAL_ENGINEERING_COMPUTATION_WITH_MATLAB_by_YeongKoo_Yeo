function Pv = vpRM(T,Tb,nu,dB,GI,Dp) 
% Estimation of vapor pressure by using the Rarey/Moller equation 
% input 
% T,Tb: temperature and normal boiling point (K) 
% nu: frequency of groups 
% dB, GI: delta B and GIij of each group 
% Dp: D prime 
% output 
% Pv: vapor pressure (bar) 
sumdBi = sum(nu.*dB); sumGI = sum(sum(GI)); Bp = 9.42208 + sumdBi + sumGI;  
Pv = exp(Bp*(T-Tb)./(T+2.65-(Tb.^1.485)/135) + Dp*log(T./Tb)); % bar 
fprintf('Vapor pressure = %g bar \n', Pv); 
end 