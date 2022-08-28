% resbutane.m: residual properties of n-butane 
eosset = {'VR', 'VDW', 'RK', 'SRK', 'PR'}; 
Tc = 425.1; Pc = 37.96; w = 0.2; T = 500; P = 50; state = 'v'; 
for i = 1:length(eosset)    
eos = eosset{i}; [Z V dH dS] = deptfun(state,eos,T,P,Tc,Pc,w);     
fprintf('The equation of state=%s: Z=%g H^R=%g S^R=%g\n', eos,Z,dH,dS); 
end 