A = [0 1;4 3;4 1]; b = [200; 1280; 960]; c = [-1 -1]; constr = '< < <';
[xopt, fopt] = simplexlp(A, b, c, constr)
