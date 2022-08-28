% 
f  =  @(x,y,z) 64*x.*y.*(1-x).^2.*z; 
I2 =  integral3(f, 0, 1, 0, 1, 0, 1) 
Id  =  triplequad(f, 0, 1, 0, 1, 0, 1) 