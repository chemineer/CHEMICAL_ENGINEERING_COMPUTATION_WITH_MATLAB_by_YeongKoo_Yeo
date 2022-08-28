c = [-2 -3]; A = [4 10;4 4;-1 0;0 -1]; b = [45 23 0 0]; intcon = [1 2];
[x, fv] = intlinprog(c, intcon, A, b)
