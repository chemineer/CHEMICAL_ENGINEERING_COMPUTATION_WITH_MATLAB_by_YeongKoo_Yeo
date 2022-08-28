% rateconst.m 
T = [313 319 323 328 333]; k = 1e-3*[0.43 1.03 1.80 3.55 7.17]; % data 
x = 1./T; y = log(k); c = polyfit(x,y,1); A = exp(c(2)), E = -c(1)*8.314 
xv = linspace(min(x),max(x),100); yv = polyval(c,xv); plot(xv,yv,x,y,'o'), 
xlabel('1/T(K^-1)'), ylabel('ln(k)') 