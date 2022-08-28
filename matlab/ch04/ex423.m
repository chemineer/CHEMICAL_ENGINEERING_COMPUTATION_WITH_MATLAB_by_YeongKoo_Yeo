% bubbleT.m : bubble temperature of a mixture 
clear all; 
% Data 
n = 4; P = 1; % n: number of components, P: total pressure (atm) 
Rg = 82.06; % gas constant (cm^3 atm/(mol K)) 
x = [0.162 0.068 0.656 0.114]; % liquid-phase mole fraction 
% Antoine equation constants 
A = [9.2033 12.2786 9.1690 9.2675]; B = [2697.55 3803.98 2731.00 2788.51]; 
C = [-48.78 -41.68 -47.11 -52.36]; 
% Virial coefficients 
Bij = zeros(n,n); Bij(1,1) = -1360.1; Bij(1,2) -657.0; Bij(1,3) = -1274.2; Bij(1,4) = -1218.8; 
Bij(2,2) = -1174.7; Bij(2,3) = -621.8; Bij(2,4) = -589.7; Bij(3,3) = -1191.9; 
Bij(3,4) = -1137.9; 
Bij(4,4) = -1086.9; 
% Data for calculation of activity coefficients using UNIFAC method 
k = 6; % number of functional groups (k) 
R = [0.9011 0.6744 0.4469 0.2195 0.5313 1.0000]; % volume vector for each functional group 
Q = [0.848 0.540 0.228 0.000 0.400 1.200]; % surface area vector for each functional group 
nu = [2 1 3 0; 4 1 1 0; 0 0 1 0; 0 0 1 0; 0 0 0 6; 0 1 0 0]; % number of functional groups 
amn = [0 0 0 0 61.13 986.5; 0 0 0 0 61.13 986.5; 0 0 0 0 61.13 986.5;...    
0 0 0 0 61.13 986.5; -11.12 -11.12 -11.12 -11.12 0 636.1;...    
156.4 156.4 156.4 156.4 89.60 0]; % interaction parameter matrix 
Terr = 1; Tcrit = 1e-6; % Initialization 
% delta(i,j) 
for i = 1:n-1, Bij(:,i) = Bij(i,:); end 
delta = zeros(n,n); 
for i = 1:n    
for j = 1:n        
delta(j,i) = 2*Bij(j,i) - Bij(i,i) - Bij(j,j);        
if j == i, delta(i,j) = 0; end    
end 
end 
for i = 1:n-1, delta(:,i) = delta(i,:); end 
% Step 1) 
PHIi = ones(1,n); % PHI_i = 1 for each component 
Tsat = B./(A - log(P)) - C; %Ti^sat 
T = sum(x.*Tsat); 
% Step 2) 
Psat = exp(A - B./(T + C)); % Pi^sat 
gamma = unifgam(k,R,Q,nu,amn,n,x,T); 
P1sat = P/(sum(x.*gamma.*Psat./PHIi/Psat(1))); % j=1 
T = B(1)/(A(1) - log(P1sat)) - C(1); 
% Step 3) 
Told = T; 
while Terr > Tcrit    
Psat = exp(A - B./(Told + C)); y = x.*gamma.*Psat./PHIi/P;    
for i = 1:n % PHI        
sumy = 0;        
for j = 1:n            
for ki = 1:n, sumy = sumy + y(j)*y(ki)*(2*delta(j,i) - delta(j,ki)); end        
end        
PHIi(i) = exp((Bij(i,i)*(P - Psat(i)) + sumy*P/2) /(Rg*Told));    
end    
gamma = unifgam(k,R,Q,nu,amn,n,x,Told);    
P1sat = P/(sum(x.*gamma.*Psat./PHIi/Psat(1))); % j=1    
T = B(1)/(A(1) - log(P1sat)) - C(1); Terr = abs(Told - T); Told = T; 
end 
T, y, Psat,PHIi, gamma 