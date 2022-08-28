%Reduction of an Iron Ore
f = @(r) 9.496e3*(1-12*r^2 + 16*r^3) - 1800; 
df = @(r) 9.496e3*(-24*r + 48*r^2); % 1st derivative of f(r) (=df/dr)
x = bisectn(f,0,0.5) % bisection method 
x = secant(f,0,0.5) % secant method 
x = newtrap(f,df,0.3) % Newton-Raphson method