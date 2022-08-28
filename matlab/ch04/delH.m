function [Q, mc] = delH(C,T1,T2)
% Calculates enthalpy change and avg. heat capacity of pure material 
% input: 
% C: coefficients of Cp relation (C=[A,B,C,D]) 
% T1,T2: temperature range (lower and upper limits of integral) 
% output: 
% Q: delta H 
R = 8.314; % J/(mol-K) 
fH = @(T) C(1) + C(2)*T + C(3)*T.^2 + C(4)*T.^(-2); %T(K) 
intCp = quadl(fH,T1,T2); Q = R*intCp; % J 
mc = intCp/(T2-T1); % (Cp)h/R 
end 