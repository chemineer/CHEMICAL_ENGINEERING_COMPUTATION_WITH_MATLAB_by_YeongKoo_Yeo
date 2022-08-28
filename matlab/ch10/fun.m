function fv = fun(x)
fv = [2*x(1) + x(2) ^2; %objective function
x(1)^2 + x(2)^2 - 8; % equality constraint: h(x) = 0
x(1); % inequality constraint: g1(x) >= 0
-x(1) + 4; % inequality constraint: g2(x) >= 0
x(2) - 1; % inequality constraint: g3(x) >= 0
-x(2) + 5]; %inequality constraint: g4(x) >= 0
end
