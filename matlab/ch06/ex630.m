% etcrack.m: ethane cracking 
% 1: C2H6, 2: C2H4, 3: C2H2, 4: C3H6; 5: C3H8, 6: C4H6, 7: CH4, 8: H2 
% Solution of differential equation system 
P = 3; T = 1073.15; Ri = 0.08314; F10 = 1300; Fs = 0.4*F10; % data 
Vspan = [0 20000]; F0 = zeros(1,8); F0(1) = F10; 
[V F] = ode45(@C2H6fun,Vspan,F0,[],P,T,Fs); n = length(V); C = []; 
for k = 1:n, Ft = sum(F(k,:)); C = [C; F(k,:)*P/(Ft*Ri*T)]; end 
% Concentration profiles 
C1 = C(:,1); C2 = C(:,2);  
plot(V,C1,V,C2,'--'), xlabel('V(liter)'), ylabel('C(mol/l)'), grid, legend('C_2H_6','C_2H_4') 
% Conversion and molal flow rate at the end of the reactor volume: 
x = (F10 - F(end,1))/F10; % conversion 
fprintf('Fractional conversion of ethane = %g\n', x); 
fprintf('Flow rate of each component at the end of the reactor volume:\n'); 
fprintf('C2H6 = %g mol/s, C2H4 = %g mol/s\n', F(end,1), F(end,2)); 
fprintf('C2H2 = %g mol/s, C3H6 = %g mol/s\n', F(end,3), F(end,4)); 
fprintf('C3H8 = %g mol/s, C4H6 = %g mol/s\n', F(end,5), F(end,6)); 
fprintf('CH4 = %g mol/s, H2 = %g mol/s\n', F(end,7), F(end,8)); 