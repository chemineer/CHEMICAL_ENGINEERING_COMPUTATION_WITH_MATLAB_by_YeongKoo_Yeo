% pipedia.m  
rho = 997.92; Q = 0.0085; mu = 0.000982; L = 250; gc = 1; Ff = 27.5; % data 
g = @(D) 0.0936*(4*rho*Q/(pi*mu*D))^(-0.2)*(L/D/gc)*(4*Q/pi/D^2)^2 - Ff; % define eqn 
D0 = 1; D = fzero(g,D0)  % calls fzero to solve the eqn. 