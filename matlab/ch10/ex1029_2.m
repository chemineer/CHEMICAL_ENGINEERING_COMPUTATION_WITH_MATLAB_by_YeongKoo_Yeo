f = @(x) (x(1)-1/2)^2*(x(1)+1)^2 + 2*(x(2)+1)^2*(x(2)-1)^2;
A = [2 4;-3 1]; b = [7; 3]; x0 = [0 0]; [x fv] = fmincon(f, x0, A, b)
