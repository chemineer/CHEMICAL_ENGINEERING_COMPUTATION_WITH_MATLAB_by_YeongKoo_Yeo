% Fugacity Coefficients in a Mixture
state = 'v'; k = zeros(2,2); eos='pr'; T = 100; P = 4.119e5; ni = [0.958 0.042]; 
Tc = [126.1 190.6]; Pc = [33.94 46.04]*1e5; w = [0.04 0.011];  
[Z,V,phi] = phimix(ni,P,T,Pc,Tc,w,k,state,eos) 