H = [2 -2;-2 4]; c = [-4; 0]; A = [2 1;1 -4;-1 0;0 -1]; b = [6 0 0 0]'; 
[x,f] = quadprog(H, c, A, b)
