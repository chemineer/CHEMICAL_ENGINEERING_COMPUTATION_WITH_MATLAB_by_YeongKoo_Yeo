% Flow in a Pipe
clear all; T=60; L=1000; D=7.981; rf=0.00015; dz=300; dP=-150;
v0 = 10; v = fzero(@vwfun, v0, [], T,L,D,rf,dz,dP), 
q = (7.481*60)*(pi*v.*(D/ 12).^2)/4 % gpm