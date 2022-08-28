% Equilibrium Conversion
P = 2; T = 340; R = 0.082; Kc = 0.1; ya0 = 1; epsilon = 1; Ca0 = ya0*P/(R*T); x0 = 0.5; 
fF = @(x) 4*Ca0*x^2 - Kc*(1-x)*(1+epsilon*x); Xe = fzero(fF,x0) 