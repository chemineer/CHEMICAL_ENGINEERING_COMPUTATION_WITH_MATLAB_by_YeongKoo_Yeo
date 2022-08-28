% wilsonpar.m: Wilson equation parameters for a binary system 
w1 = 78; mw1 = 46.07; mw2 = 114; Ta = 77; P = 760; 
A1 = 8.04494; B1 = 1554.3; C1 = 222.65; A2 = 6.92374; B2 = 1355.126; C2 = 209.517; 
x1 = (w1/mw1)/(w1/mw1 + (100-w1)/mw2); x2 = 1 - x1; % mole fraction 
P1 = 10^(A1 - B1/(Ta + C1)); P2 = 10^(A2 - B2/(Ta + C2)); % vapor pressure 
gam1 = P/P1; gam2 = P/P2; % activity coefficients 
% solve the nonlinear system (g(1)=g12, g(2)=g21) 
g0 = [0.1 0.1]; g = fsolve(@wilact,g0,[],x1,x2,gam1,gam2); 
g12 = g(1), g21 = g(2)  

function f = wilact(g,x1,x2,gam1,gam2) 
t1 = x1 + g(1)*x2; t2 = x2 + g(2)*x1; 
f = [log(gam1) + log(t1) - (g(1)*t2 - g(2)*t1)*x2/(t1*t2);    
log(gam2) + log(t2) + (g(1)*t2 - g(2)*t1)*x1/(t1*t2)]; 
end