function hum = humidest(Hr,Td) 
% Hr: relative humidity 
% Td: dry bulb temperature  

% critical properties 
Tc = 647.096; % critical temperature (K) 
Pc = 22064000; % critical pressure (Pa) 
% constants 
Rw = 461.512244565; 
a = [-7.85951783 1.84408259 -11.7866497 22.6807411 -15.9618719 1.80122502]; 
% temperature parameters 
Td = Td+273.15; theta = Td/Tc; tau = 1 - theta; 
% saturated pressure 
tw = (Tc/Td)*(a(1)*tau + a(2)*tau^1.5 + a(3)*tau^3+a(4)*tau^3.5 + a(5)*tau^4 + a(6)*tau^7.5); 
Ps = Pc*exp(tw); 
% results 
hum.pres =(Hr/100)*Ps; % actual vapor pressure 
hum.Ha = hum.pres*1000/((Td)*Rw); % absolute humidity 
fprintf('Actual vapor pressure: %g Pa\n', hum.pres); 
fprintf('Absolute humidity = %g g/m^3\n', hum.Ha); 