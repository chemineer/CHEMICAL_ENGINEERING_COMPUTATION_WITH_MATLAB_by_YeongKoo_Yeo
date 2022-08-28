function dF = C2H6fun(V,F,P,T,Fs) 
% 1: C2H6, 2: C2H4, 3: C2H2, 4: C3H6; 5: C3H8, 6: C4H6, 7: CH4, 8: H2 
Re = 1.987; Ri = 0.08314; % gas constant 
k1 = 4.65e13*exp(-65210/(Re*T)); k2 = 3.85e11*exp(-65210/(Re*T)); 
k3 = 9.81e8*exp(-36920/(Re*T));  k4 = 1.03e12*exp(-41260/(Re*T)); 
k5 = 7.08e13*exp(-60430/(Re*T)); 
Ft = sum(F) + Fs; C = F*P/(Ft*Ri*T);  
r1 = k1*C(1); r2 = k2*C(1)^2; r3 = k3*C(4); r4 = k4*C(2)*C(3); r5 = k5*C(1)*C(2); 
dF = [-r1 - 2*r2 - r5; r1 - r4 - r5; r3 - r4; -r3 + r5; r2; r4; r2 + r3 + r5; r1]; 
end 