% filtintp.m: cubic spline interpolation by the built-in spline function 
x  =  [0 9.7 15.6 21.3 31.7 35.2 38.4 42.9]; 
y  =  [0 0.28 0.524 0.998 1.695 2.306 2.781 3.205]; 
xi  =  min(x):0.1:max(x); yi  =  spline(x,y,xi); 
plot(xi,yi,x,y,'o'), xlabel('Flow rate(liter/sec)'), ylabel('Pressure drop(kPa)') 
legend('Cubic spline interpolation','Experimental data','Location','Best')