%Nonlinear Equation System
f =  @(x) [cos(x(1)) + x(2)^2 + log(x(3)) - 8; 4*x(1) + 3^x(2) - x(3)^3 + 2;x(1) + x(2) + x(3) - 6]; 
x0  =  [1 1 1]'; z  =  newtrapmv(f,x0)