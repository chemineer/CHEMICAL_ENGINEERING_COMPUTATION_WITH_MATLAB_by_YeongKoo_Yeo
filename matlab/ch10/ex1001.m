%Fibonacci Search
a = -8; b = 8; n = 20; fobj = @(x) 2*x^2*sin(1.5*x) - 4*x^2 + 3*x-1; 
[x, f, fint] = fbnopt(fobj, a, b, n)