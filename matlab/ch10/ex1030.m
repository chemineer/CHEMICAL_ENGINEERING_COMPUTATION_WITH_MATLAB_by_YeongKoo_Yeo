f = [-50 -98 -25 -43]'; A = [6 12 3 8;4 30 2 1;-1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 -1]; b = [1150 750 0 0 0 0]';
[x, fv] = linprog(f, A, b)
