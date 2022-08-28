% mmmdl.m: Michaelis-Menten equation 
S  =  [1.2 1.6 3.2 4.3 5.8 7.6 8.8]; r  =  [0.06 0.12 0.24 0.27 0.33 0.34 0.34]; 
p  =  polyfit(1./S,1./r,1); % p(1) = p1, p(2) = p0 
a  =  1/p(2), b  =  p(1)*a 
Sv =  S(1):0.1:S(end); rv =  a*Sv./(b + Sv); % Reaction rates by Michaelis-Menten model 
plot(Sv,rv,S,r,'o'),xlabel('[S]'),ylabel('r'),legend('Michaelis-Menten model','Experimental data')  