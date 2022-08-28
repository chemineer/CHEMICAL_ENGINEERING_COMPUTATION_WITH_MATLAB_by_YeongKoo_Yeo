function metT 
y0 = [100 0]; metin = bvpinit(linspace(0,1,20),y0); 
% divide the interval [0 1] into 20 subintervals 
y  =  bvp4c(@metfun, @metbc, metin); 
xv  =  linspace(0,1); yv  =  deval(y,xv); plot(xv,yv(1,:)), xlabel('x'), ylabel('T(degC)'), grid 
end 
% Define differential equations 
function dy  =  metfun(x,y) 
h  =  50; D  =  0.04; k  =  390; Ta  =  25; 
dy  =  [y(2); 4*h*(y(1)-Ta)/(D*k)]; 
end 
% Define boundary conditions 
function mb  =  metbc(ya,yb) 
mb  =  [ya(1)-100; yb(1)]; 
end 