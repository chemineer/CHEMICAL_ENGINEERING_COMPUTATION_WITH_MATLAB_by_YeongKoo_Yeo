A = [2 5; 1 1; 3 1]; b = [20; 6; 9]; c = [2 1]; constr = '> > >';
[xopt, fopt] = simplexlp(A, b, c, constr)
