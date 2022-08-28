% Double integral 
fxy  =  @(x,y) 3*x.^y+x-1.2*x.^2 -3*y.^2+25; 
I2  =  integral2(fxy, 0, 7, 0, 5) 
Id  =  dblquad(fxy, 0, 7, 0, 5)