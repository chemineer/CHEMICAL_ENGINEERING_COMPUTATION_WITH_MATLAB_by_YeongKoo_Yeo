% btbubbleT.m: bubble point temperature of beznene-toluene mixture 
% 1: benzene, 2: toluene 
A  =  [15.90085, 16.01066]; B  =  [2788.507, 3094.543]; C  =  [220.790, 219.377]; 
P  =  760; % total pressure 
x  =  [0.4, 0.6]; % mole fraction 
f  =  @(t) x(1)*exp(A(1)-B(1)/(t+C(1))) + x(2)*exp(A(2)-B(2)/(t+C(2))) - P; 
t  =  fzero(f,50); % initial guess  =  50 deg.C 
fprintf('The bubble point temperature  =  %g deg.C\n', t);