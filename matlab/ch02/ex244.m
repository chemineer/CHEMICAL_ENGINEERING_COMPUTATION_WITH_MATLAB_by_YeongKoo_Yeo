% Use of integral, quad, and quadl 
format long 
q  =  0.2; r  =  0.7; s  =  4; 
hfun  =  @(x) 1./((x-q).^2+0.01)+1./((x-r).^2+0.04)-s; 
integral(hfun,0,1)
quad(hfun,0,1)
quadl(hfun,0,1)