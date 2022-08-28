% dfgrad.m 
f  =  @(x) 0.3+20*x-180*x.^2+650*x.^3 -880*x.^4+360*x.^5; 
df  =  @(x) 20 -360*x+1950*x.^2 -3520*x.^3+1800*x.^4; 
x  =  0:0.1:1.0; y  =  f(x); 
dr  =  gradient(y, 0.1), dy  =  df(x) % dr: estimation by gradient, dy: exact solution 