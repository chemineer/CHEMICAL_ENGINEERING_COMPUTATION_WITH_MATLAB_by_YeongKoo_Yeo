% enthmix.m 
state = 'L'; nx = [0.419 0.3783 0.2027]; P = 6.8947e5; T = 158; eos = 'rk';  
Pc = [45.99 48.72 42.48]*1e5; Tc = [190.6 305.3 369.8]; w = [0.012 0.1 0.152]; 
k = zeros(3,3); 
Afi = [8.245223 0.3806333e-2 0.8864745e-5 -0.7461153e-8 0.182296e-11;    
11.51606 0.140309e-1 0.854034e-5 -0.1106078e-7 0.31622e-11;    
15.58683 0.2504953e-1 0.1404258e-4 -0.352626e-7 0.1864467e-10]; 
state = upper(state); if state == 'L', x = nx; else x = ny; end  
[Z H] = khmix(x,P,T,state,eos,Pc,Tc,w,k,Afi); 
fprintf('Equation of state: %s, State: %s',upper(eos),upper(state)); 
fprintf('\nCompressibility factor = %g, Enthalpy = %g (J/mol)\n',Z,H);