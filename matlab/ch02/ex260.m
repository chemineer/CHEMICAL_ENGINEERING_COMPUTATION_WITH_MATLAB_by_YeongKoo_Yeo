% sstem.m: steady-state temperature distribution 
% u_xx + u_yy  =  0 
R  =  [0 4 0 4]; f  =  @(x,y) 0; g  =  @(x,y) 0; 
qx0  =  @(y) exp(y) - cos(y); qxf  =  @(y) exp(y)*cos(4) - exp(4)*cos(y); 
qy0  =  @(x) cos(x) - exp(x); qyf  =  @(x) exp(4)*cos(x) - exp(x)*cos(4); 
m  =  40; n  =  40; crit  =  1e-6; kmax  =  1000; 
[u,x,y]  =  helpde(f,g,qx0,qxf,qy0,qyf,R,m,n,crit,kmax); 
mesh(x,y,u), xlabel('x'), ylabel('y'), zlabel('u(x,y)'), colormap(gray), 
axis([0 4 0 4 -100 100]) 