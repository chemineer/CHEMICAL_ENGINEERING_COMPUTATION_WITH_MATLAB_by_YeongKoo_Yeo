% rmin.m
zf = [0.0137 0.5113 0.0411 0.0171 0.0262 0.0446 0.3106 0.0354];
relvol = [2.43 1.93 1.00 0.765 0.362 0.164 0.0720 0.0362];
d = [12 442 13 0 0 0 0 0]; q = 0.87; x0 = [1.5 0.8 3 500 1.5];
x = fsolve(@debutanizer, x0, [], zf, q, relvol, d)
function f = debutanizer(x, zf, q, relvol, d)
% zf: feed composition, q: feed quality, relvol: relative volatility, d: component rate in D
% x(1): theta1, x(2): theta2, x(3): d(4), x(4): D, x(5): Rmin
d(4) = x(3);
f = [sum(d) - x(4);
sum(relvol.*zf./(relvol - x(1))) - (1-q);
sum(relvol.*zf./(relvol - x(2))) - (1-q);
sum(relvol.*d./(relvol - x(1))) - x(4)*(1 + x(5));
sum(relvol.*d./(relvol - x(2))) - x(4)*(1 + x(5))];
end
