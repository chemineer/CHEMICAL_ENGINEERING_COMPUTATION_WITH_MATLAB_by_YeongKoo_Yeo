%Redlich-Kwong equation
T = 350; P = 9.4573; Tc = 425.1; Pc = 37.96; w = 0.2; state = 'L'; eos = 'rk';  
[Z V] = cubicEOSZ(state,eos,T,P,Tc,Pc,w)