function [f g] = grgfun(x)
x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);
x5 = x(5); x6 = x(6); x7 = x(7); % slack variables
f = x1^2 + x2^2 + 2.3*x3^2 - 1.2*x4^2 - 4*x1 - 6*x2 - 20*x3 + 6*x4 + 100;
g(1) = x1^2 + x2^2 + x3^2 + x4^2 + x1 - x2 + x3 - x4 + x5 - 7;
g(2) = x1^2 + 2*x2^2 + x3^2 + 2*x4^2 - x1 - x4 + x6 - 11;
g(3) = 2*x1^2 + x2^2 + x3^2 + 2*x1 - x2 - x4 + x7 - 6;
end
