% flashbubP.m : flash calculations for multicomponent system (Bubble P) 
clear all; 
% Data 
n = 4; P = 1; T = 334.15; % n: number of components, P: atm, T: K 
Rg = 82.06; % gas constant (cm^3 atm/(mol K)) 
z = [0.25 0.40 0.20 0.15]; % feed stream mole fraction 
% Antoine equation constants 
A = [9.2033 12.2786 9.1690 9.2675]; B = [2697.55 3803.98 2731.00 2788.51]; 
C = [-48.78 -41.68 -47.11 -52.36]; 
% Virial coefficients 
Bij = zeros(n,n); 
Bij(1,1) = -1360.1; Bij(1,2) -657.0; Bij(1,3) = -1274.2; Bij(1,4) = -1218.8; 
Bij(2,2) = -1174.7; Bij(2,3) = -621.8; Bij(2,4) = -589.7; 
Bij(3,3) = -1191.9; Bij(3,4) = -1137.9; Bij(4,4) = -1086.9; 
% Data for calculation of activity coefficients using UNIFAC method 
k = 6; % number of functional groups (k) 
R = [0.9011 0.6744 0.4469 0.2195 0.5313 1.0000]; % volume vector for each functional group 
Q = [0.848 0.540 0.228 0.000 0.400 1.200]; % surface area vector for each functional group 
nu = [2 1 3 0; 4 1 1 0; 0 0 1 0; 0 0 1 0; 0 0 0 6; 0 1 0 0]; % number of functional groups 
amn = [0 0 0 0 61.13 986.5; 0 0 0 0 61.13 986.5; 0 0 0 0 61.13 986.5;...    
0 0 0 0 61.13 986.5; -11.12 -11.12 -11.12 -11.12 0 636.1;...    
156.4 156.4 156.4 156.4 89.60 0]; % interaction parameter matrix 
Aerr = 1; Acrit = 1e-6; % Initialization 
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
x = z; % Bubble P 
Psat = exp(A - B./(T + C)); % Pi^sat 
gamma = unifgam(k,R,Q,nu,amn,n,x,T); 
PHIi = ones(1,n); % initial guess 
y = x.*gamma.*Psat./PHIi/P; 
for i = 1:n % PHI    
sumy = 0;    
for j = 1:n        
for m = 1:n, sumy = sumy + y(j)*y(m)*(2*delta(j,i) - delta(j,m)); end    
end    
PHIi(i) = exp((Bij(i,i)*(P - Psat(i)) + sumy*P/2) /(Rg*T)); 
end 
% Step 2) 
K = gamma.*Psat./(PHIi*P); alpha0 = 0.5; alphaold = alpha0; 
% Step 3) 
while Aerr > Acrit    
falpha = @(alpha) sum((1-K).*z./(1 + alpha*(K - 1)));    
alpha = fzero(falpha, alphaold); x = z./(1 + alpha*(K - 1)); y = K.*x;    
gamma = unifgam(k,R,Q,nu,amn,n,x,T);    
for i = 1:n % PHI        
sumy = 0;        
for j = 1:n            
for m = 1:n, sumy = sumy + y(j)*y(m)*(2*delta(j,i) - delta(j,m)); end        
end        
PHIi(i) = exp((Bij(i,i)*(P - Psat(i)) + sumy*P/2) /(Rg*T));    
end    
K = gamma.*Psat./(PHIi*P); Aerr = abs(alpha - alphaold); alphaold = alpha; 
end 
alpha, x, y, K, gamma 