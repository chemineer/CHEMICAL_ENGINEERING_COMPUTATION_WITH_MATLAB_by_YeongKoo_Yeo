% Water-Isobutanol Equilibrium Calculations
x110 = 0; x220 = 0; x210 = 1; x120 = 1; T0 = 100; beta0 = 0.8; 
z = [0.2 0.8]; P = 760; t0 = [x110 x210 x120 x220 T0 beta0]; 
[x,fval] = fsolve(@beq, t0, [], z, P); x  