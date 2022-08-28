% qdbvp4c.m 
% 'bvp4c' undefined
initsol  =  bvpinit(linspace(1,2,10), [1 0]); res  =  bvp4c(@bcprob, @bcval, initsol); 
x  =  linspace(1,2); y  =  deval(res,x); plot(x,y(1,:)), grid, xlabel('x'), 
ylabel('y') 
% define differential equations 
function deset  =  bcprob(x,z) 
deset  =  [z(2); 6*z(1)./x.^2]; 
end 
% boundary conditions 
function bcset  =  bcval(za,zb) 
bcset  =  [za(1)-1; zb(1)-1]; 
end 