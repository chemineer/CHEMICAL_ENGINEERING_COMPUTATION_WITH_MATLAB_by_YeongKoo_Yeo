% wsprop.m: calculate properties (H,S,U,V) of H2O (water and steam) using H2Oprop.m 
P = input('Pressure (kPa) = '); 
T = input('Temperature (deg.C) = '); 
w = H2Oprop(P,T); 
fprintf('Specific enthalpy = %g kJ/kg\n', w.H) 
fprintf('Specific entropy = %g kJ/(kg-K)\n', w.S) 
fprintf('Specific internal energy = %g kJ/kg\n', w.U) 
fprintf('Specific volume = %g cm^3/g\n', w.V)