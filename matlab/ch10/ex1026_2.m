f = @(x) 2*x^2*sin(x) + exp(-x); lb = -4; ub = 0; 
[x,fv] = fminbnd(f, lb ,ub)
