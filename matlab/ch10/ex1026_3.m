f = @(x) 2*x^2*sin(x) + exp(-x); lb = -4; ub = 0; x0 = -1; x = lsqnonlin(f, x0, lb, ub)
f = @(x) (2*x^2*sin(x) + exp(-x))^2; lb = -4; ub = 0; x0 = -1; [x,fv] = fminbnd(f, lb ,ub)
