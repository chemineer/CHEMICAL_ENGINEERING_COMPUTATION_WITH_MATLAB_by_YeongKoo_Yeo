% bvp4cex.m 

zinit  =  bvpinit(linspace(0,4,10),[1 0]); sol  =  bvp4c(@zode,@zobc,zinit); 
x  =  linspace(0,4,100); y  =  deval(sol,x); plot(x,y(1,:)), grid, xlabel('x'), ylabel('y') 
% define differential equations 
function dz  =  zode(x,z) 
dz  =  [z(2); -abs(z(1))]; 
end 
% boundary conditions 
function res  =  zobc(za,zb) 
res  =  [za(1); zb(1)+2]; 
end 