% bte1dintp.m: 1D interpolation for B-T equilibrium data 
x  =  [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0]; 
y  =  [0.00 0.21 0.38 0.51 0.62 0.71 0.79 0.86 0.91 0.96 1.00]; 
lv  =  interp1(x,y,0.45,'linear'); % Linear interpolation 
pv  =  interp1(x,y,0.45,'pchip'); % Piecewise cubic Hermite interpolation 
sv  =  interp1(x,y,0.45,'spline'); % Piecewise cubic spline 
nv  =  interp1(x,y,0.45,'nearest'); % Nearest neighbor interpolation 
fprintf('linear: %6.4f\npchip: %6.4f\nspline: %6.4f\nnearest: %6.4f \n',lv,pv,sv,nv); 