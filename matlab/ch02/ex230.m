% enth2nd.m: interpolation of enthalpy for saturated steam 
format long; 
T  =  [283.15 303.15 323.15 363.15 393.15 413.15]; % T(K) 
H  =  [2519.9 2556.4 2592.2 2660.1 2706.0 2733.1]; % enthalpy (kJ/kg) 
p  =  polyfit(T, H, 2) 
Tv  =  280:0.1:415; Hv  =  polyval(p,Tv); % enthalpy values by the 2nd polynomial 
plot(Tv,Hv,T,H,'o'), xlabel('T(K)'), ylabel('H(kJ/kg)') 
legend('2nd-order interpolation','Steam table','location','best')  
f  =  polyval(p, 350.15)