% Numerical Integration
f  =  @(x) 1/(1+x^2); a  =  0; b  =  1.5; n  =  10;
z  =  trapzoid(f,a,b,n) % trapezoidal rule
z  =  simps(f,a,b,n) % Simpson 1/3 rule