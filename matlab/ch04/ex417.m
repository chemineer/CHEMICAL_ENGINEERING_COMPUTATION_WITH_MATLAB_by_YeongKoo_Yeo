% Estimation of Pressure by Raoult’s Law
v = 0.65; T = 60; z = 0.6; A = [13.8183 13.8587]; B = [2477.07 2991.32]; 
C = [233.21 216.64]; 
x0 = [0.1 0.6 50]; x = fsolve(@binflash, x0, [], v, T, z, A, B, C) 