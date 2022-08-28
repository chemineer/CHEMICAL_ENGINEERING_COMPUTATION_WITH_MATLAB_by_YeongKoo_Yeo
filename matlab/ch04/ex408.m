% Enthalpy and Entropy Departures of Propane Gas
Tc = 369.8; Pc = 42.49; w = 0.152; T = 378.15; P = 5; eos = 'pr'; state = 'v'; 
T1 = 378.15; P1 = 5; T2 = 463.15; P2 = 25; 
[Z1 V1 dH1 dS1] = deptfun(state,eos,T1,P1,Tc,Pc,w) % state 1 
[Z2 V2 dH2 dS2] = deptfun(state,eos,T2,P2,Tc,Pc,w) % state 2
A = -4.224; B = 0.3063; C = -1.586e-4; D = 3.215e-8; Tc = 369.8;  
Pc = 42.49; w = 0.152;  
eos = 'pr'; state = 'v'; T1 = 378.15; P1 = 5; T2 = 463.15; P2 = 25; 
[dH dS] = delHS(state,eos,T1,P1,T2,P2, A,B,C,D,Tc,Pc,w)