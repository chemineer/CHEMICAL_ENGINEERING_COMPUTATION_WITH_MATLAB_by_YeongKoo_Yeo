% Equilibrium of Water-Gas Shift Reaction
f = @(x) 148.4 - x^2/(1 - x)^2; x0 = 0.5; x = fzero(f,x0)  
% initial guess = 0.5. Uses the built-in solver fzero. 