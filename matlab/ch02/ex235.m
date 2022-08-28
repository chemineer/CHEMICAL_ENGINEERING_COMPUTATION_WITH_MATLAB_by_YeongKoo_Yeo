% Two-Dimensional Interpolation

x  =  [2 10]; y  =  [1 8]; z  =  [80 78; 75 90]; % data 
zi  =  interp2(x, y, z, 6.4, 5.2, 'spline') 