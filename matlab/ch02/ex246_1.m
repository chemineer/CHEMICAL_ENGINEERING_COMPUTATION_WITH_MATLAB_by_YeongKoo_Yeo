% Double integral 
f  =  @(x,y) 1./(1-x.*y); 
I2 =  integral2(f, 0, 0.8, 0, 0.8) 
I2 =  dblquad(f, 0, 0.8, 0, 0.8)