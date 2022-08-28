%Friction Factor Using the Colebrook Equation 
eD  =  1.3e-4; Nre  =  6.5e4; f0  =  0.1; % data (f0: initial guess) 
Cf =  @(f) 1/sqrt(f) + 0.86*log(eD/3.7 + 2.51/Nre/sqrt(f)); % define equation
f  =  fzero(Cf,f0)