% regorder.m 
t = [0 50 100 150 200 250 300]; Ca = [0.05 0.038 0.0306 0.0256 0.0222 0.0195 0.0174]; % data 
Ca0 = 0.05; x = fsolve(@frx,[0.2 2],[],t,Ca0,Ca) % initial value: k=0.2, alpha=2 

function f = frx(x,t,Ca0,Ca) 
% x(1)=k’, x(2)=alpha 
n = length(t);  
for i = 1:n, f(i,1) = x(1)*(1-x(2))*t(i) - Ca0^(1-x(2)) + Ca(i)^(1-x(2)); end 
end 