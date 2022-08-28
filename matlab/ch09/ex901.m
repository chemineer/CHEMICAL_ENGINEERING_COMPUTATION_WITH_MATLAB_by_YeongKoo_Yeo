%'pkg load symbolic' required
syms t a b; 
f = 1 + t + t^2 + sin(a*t) - t*cos(b*t); Lf = laplace(f) 