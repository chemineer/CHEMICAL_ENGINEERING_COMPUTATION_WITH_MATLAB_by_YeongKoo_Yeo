% bubbp.m: calculation of total P and T 
clear all; 
T = 75; A1 = 14.2724; B1 = 2945.47; C1 = 224; A2 = 14.2043; B2 = 2972.64; C2 = 209; 
fp = @(x) exp(A2-B2/(T+C2)) + x*(exp(A1-B1/(T+C1)) - exp(A2-B2/(T+C2))); 
% total pressure 
x1 = 0:0.01:1; y1 = []; P = []; 
for k = 1:length(x1) 
Ps = feval(fp,x1(k)); fy = @(x) x*(exp(A1-B1/(T+C1)))/Ps; 
y = feval(fy,x1(k)); y1 = [y1 y]; P = [P Ps]; 
end 
plot(x1,P,y1,P,'.-'), xlabel('x_1,y_1'), ylabel('P') 
legend('P-x_1','P-y_1','location','best') 