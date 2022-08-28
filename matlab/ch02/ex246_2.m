% 
f  =  @(x,y) x.^2.*y; 
I2 =  integral2(f, 1, 2, @(x) x.^2, @(x) x.^4) 
Id  =  quad2d(f, 1, 2, @(x) x.^2, @(x) x.^4) 