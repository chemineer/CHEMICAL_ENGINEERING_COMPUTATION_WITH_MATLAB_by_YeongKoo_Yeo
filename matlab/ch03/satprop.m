% satprop.m: properties of saturated H2O (vapor and liquid)  

% Basic syntax: w = satH2Oprop(dtype, phH2O, dvalue) 
% Inputs: 
%       dtype: data type( 'P' = pressure, 'T' = temperature) 
%       phH2O: phase of H2O ('V' = vapor, 'L' = liquid) 
%       dvalue: value of pressure(kPa) or temperature(deg.C) 
% Output: 
%       w- structure of property values (T,P,H,S,U, or V) of saturated H2O  
%       w.P: saturation pressure (kPa) 
%       w.T: saturation temperature (deg.C) 
%       w.H: specific enthalpy (kJ/kg) 
%       w.S: specific entropy (kJ/(kg-K)) 
%       w.U: specific internal energy (kJ/kg) 
% w.V: specific volume (cm^3/g) 

dtype = input('Data type (p(Pressure, kPa) or t(temperature, deg.C)) = '); 
dvalue = input('Data value = '); 
phH2O = input('Phase (V: vapor, L: liquid) = '); 
w = satH2Oprop(dtype, phH2O, dvalue); 
fprintf('\nSaturation pressure = %g kPa\n', w.P) 
fprintf('Saturation temperature = %g deg.C\n', w.T) 
fprintf('Specific enthalpy = %g kJ/kg\n', w.H) 
fprintf('Specific entropy = %g kJ/(kg-K)\n', w.S) 
fprintf('Specific internal energy = %g kJ/kg\n', w.U) 
fprintf('Specific volume = %g cm^3/g\n', w.V) 

% x = satH2Oprop('P','L',10540)
% x = satH2Oprop('T','V',198)
