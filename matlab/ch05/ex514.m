% powfun.m 
rho = 62.4; mu = 6.905e-4; D = 1/3; Q = 0.2; K = 0.45+3*0.5+1; gc = 32.174; 
L = 5+300+100+120+20; eD = 5e-5; 
v = 4*Q/(pi*D^2); m = Q*rho; Nre = D*v*rho/mu; % Reynolds 수 
if Nre < 2100    
f = 16./Nre; 
else % Shacham eqn    
den = 16*(log10(eD/3.7 - 5.02*log10(eD/3.7+14.5./Nre)./Nre)).^2; f = 1./den; 
end 
dH = (4*f*L/D + K)*v^2/(2*gc); Ws = m*(dH + (105-20)); 
v,f,Ws 