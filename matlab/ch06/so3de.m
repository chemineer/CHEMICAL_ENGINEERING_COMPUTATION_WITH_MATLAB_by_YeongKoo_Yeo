function dz = so3de(w,z,sa) 
% The reaction model is based on the following references: 
% [1] Harrer, T. S., Kirk Othmer encyclopedia of chemical technology,  
%     2nd ed., Vol.19, Wiley-Interscience, New York, NY, p.470, 1969. 
% [2] Mariano M. Martin (Editor), Introduction to software for chemical  
%     engineers, CRC Press, Taylor & Francis Group, Boca Raton, FL, p.136- 145, 2015. 
% Retrieve data 
Ta = sa.Ta; T0 = sa.T0; Pt0 = sa.Pt0; rho0 = sa.rho0; rhob = sa.rhob;  
ya0 = sa.ya0; yb0 = sa.yb0; yc0 = sa.yc0; Ft0 = sa.Ft0; G = sa.G;  
epn = sa.epn; phi = sa.phi; mu = sa.mu; D = sa.D; Dp = sa.Dp; U = sa.U; 
% Assign variables 
x = z(1); T = z(2); P = z(3); 
% Rate and equilibrium constants  
k = 9.8692e-3 * exp(-1.76008e5/(1.8*(T - 273.15) + 491.67) -...    
110.1*log((1.8*(T-273.15)+491.67)) + 912.8); 
Kp = 3.1415e-3 * exp(42311/(1.987*(1.8*(T-273.17) + 491.67)) - 11.24);  
% parameters (a: SO2, b: O2, c: N2, d: SO3) 
phib = yb0/ya0; Pa0 = Pt0*ya0; Fa0 = Ft0*ya0; Ac = pi*D^2/4; 
Cpsum = 300.85 - 0.0402*T + 1.8e-4*T^2 - 9.071e-8*T^3; 
dCp = -21.535 + 0.0789*T - 7.112e-5*T^2 + 2.447e-8*T^3; 
dHr = - 98480 - 21.535*(T-298) + 0.0395*(T^2-298^2) -...      
2.371e-5*(T^3-298^3) + 6.11675e-9*(T^4-298^4); 
if x < 0.05    
xs = 0.05;    
r = k*sqrt((1-xs)/xs)*(Pa0*((1.1-xs/2)/(1+epn*xs))*P/Pt0 - (xs/(Kp*(1- xs)))^2); 
else    
r = k*sqrt((1-x)/x)*(Pa0*((1.1-x/2)/(1+epn*x))*P/Pt0 - (x/(Kp*(1-x)))^2); 
end     
dz(1,1) = r/Fa0; 
dz(2,1) = (4*U*(Ta-T)/(rhob*D) + r*(-dHr))/(Fa0*(Cpsum + x*dCp)); 
dz(3,1) = -G*(1-phi)*(1+epn*x)*Pt0*T*(150*(1-phi)*mu/Dp + 1.752*G) / (P*T0*...          
rhob*Ac*rho0*Dp*phi^3); 
end 