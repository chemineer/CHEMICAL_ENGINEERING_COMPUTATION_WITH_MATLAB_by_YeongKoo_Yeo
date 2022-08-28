%A System of Nonlinear Equations
fun =  @(x) [sin(x(1))+x(2)^2+log(x(3))-7; 3*x(1)+2*x(2)-x(3)^3+1; x(1)+x(2)+x(3)-5];
x0  =  [0 2 2]'; x  =  fsolve(fun,x0) % use the built-in solver fsolve 
