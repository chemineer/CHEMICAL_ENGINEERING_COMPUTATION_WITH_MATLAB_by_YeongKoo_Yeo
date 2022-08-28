% qrmT.m: calculates properties at various T 
clear all; T=[40 60 100]; nT = length(T); L=1000; D=7.981; rf=0.00015; dz=300; 
dP=-150; v0 = 10;  
rhor = []; mur = []; qr = []; 
for i = 1:nT    
v = fzero(@vwfun, v0, [], T(i),L,D,rf,dz,dP); % calls fzero to solve the eqn    
q = (7.481*60)*(pi*v.*(D/12).^2)/4;     
rho = 62.122+0.0122*T(i)-(1.54e-4)*T(i).^2+(2.65e-7)*T(i).^3-(2.24e-10)*T(i).^4;    
mu = exp(-11.0318 + 1057.51./(T(i)+214.624));    
rhor = [rhor rho]; mur = [mur mu]; qr = [qr q];  
end 
qr, rhor, mur 