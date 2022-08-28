% isotplr.m 
v0 = 0.5; k = 0.3; Vf = 2.4; Vi = [0 Vf]; C0 = [2 0 2]; % data 
dC = @(V,C) [-2*k*C(1)^2/v0; k*C(1)^2/v0; 0]; [V C] = ode45(dC,Vi,C0); 
plot(V,C(:,1),V,C(:,2),'.-',V,C(:,3),'--'), legend('C_A','C_B','C_C') 
xlabel('V(length, m)'), ylabel('Concentration(kmol/m^3)') 