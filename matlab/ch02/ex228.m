% Determination of Interpolating 2nd-Order Polynomial
format long
T  =  [313.15 363.15 413.15]'; A  =  [T.^2 T ones(3,1)]; b  =  [8.2583 7.4799 6.9284]';
x  =  A\b
s  =  x(1)*373.15^2 + x(2)*373.15+x(3)