%SRK Equation of State 
% srkeqn.m: solves the SRK equation
R  =  0.082054; Tc  =  419.6; Pc  =  396.743; T  =  415; P  =  207.2538; % data 
a  =  0.42748*R^2*Tc^2.5/Pc; b  =  0.08664*R*Tc/Pc; % coefficients a and b 
f  =  @(V) R*T./(V-b) - a./(V.*(V+b)*sqrt(T)) - P; % define the equation 
V0  =  0.1; V  =  fzero(f,V0) % use the built-in solver fzero (V0: initial guess)  