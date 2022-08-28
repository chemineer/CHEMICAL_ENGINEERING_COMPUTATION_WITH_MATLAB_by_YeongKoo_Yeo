% fxdiff.m: differentiation using diff function 
f  =  @(x) 0.3+20*x-180*x.^2+650*x.^3 - 880*x.^4+360*x.^5; % define function 
x  =  0:0.1:1.0; n  =  length(x); % range of x 
y  =  f(x); dr  =  diff(y)./diff(x) % use diff to find f'(x) numerically 
xm  =  (x(1:n-1) + x(2:n))/2; % the midpoint of each interval 
xp  =  0:0.01:1; yp  =  20 -360*xp+1950*xp.^2 -3520*xp.^3+1800*xp.^4; % exact differentiation 
plot(xp,yp,xm,dr,'o'), xlabel('x'), ylabel('df(x)/dx') 
legend('Exact differentiation','Numerical differentiation','location','best')