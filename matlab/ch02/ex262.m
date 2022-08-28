% temprof.m 
m  =  0; x  =  linspace(0,1,20); t  =  linspace(0,21600,54); 
sol  =  pdepe(m,@pdeTde,@pdeTic,@pdeTbc,x,t); 
T = sol(:,:,1); surf(x,t,T), xlabel('x'), ylabel('t(sec)'), zlabel('T(deg C)') 
function [g,f,r]  =  pdeTde(x,t,u,DuDx) % define g, f, and r 
alpha  =  4.8e-7; g  =  1/alpha; f  =  DuDx; r  =  0; 
end 
function u0  =  pdeTic(x) % define initial conditions
u0  =  90; 
end 
function [pl,ql,pr,qr]  =  pdeTbc(xl,ul,xr,ur,t) % define boundary conditions 
pl  =  ul-15; ql  =  0; pr  =  ur-15; qr  =  0; 
end 