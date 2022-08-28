% dewpt.m: calculation of dew point 
clear all;  
A1 = 14.2724; B1 = 2945.47; C1 = 224; A2 = 14.2043; B2 = 2972.64; C2 = 209;  
T0 = 60; P = 70; x = 0:0.01:1; T = []; y = []; 
for k = 1:length(x)    
x1 = x(k);     
fp = @(T) exp(A2-B2/(T+C2)) + x1*(exp(A1-B1/(T+C1)) - exp(A2-B2/(T+C2)))-P;    
T1 = fzero(fp,T0); y1 = x1*exp(A1-B1/(T1+C1))/P; T = [T T1]; y = [y y1]; 
end 
plot(x,T,y,T,'.-'), xlabel('x_1, y_1'), ylabel('T(C)') 
legend('T-x_1','T-y_1','location','best') 