% binVLE.m: VLE calculation for binary systems using modified Raoult's law 
% Antoine vapor pressure equation parameters 
A = [16.59158 14.25326]; B = [3643.31 2665.54]; C = [33.424 53.424]; 
% (1) Bubble P calculation 
T = 318.15; x1 = 0.25; x = [x1 1-x1]; % data 
gamma = gamma12(x1,T); Psat = vp12(A,B,C,T); 
P = sum(x.*gamma.*Psat); y = (x.*gamma.*Psat)/P; 
fprintf('(1) Bubble P: P = %g, y1 = %g, y2 = %g\n',P,y(1),y(2)); 
fprintf('gamma1 = %g, gamma2 = %g\n',gamma(1),gamma(2)); 
% (2) Dew P calculation 
T = 318.15; y1 = 0.60; y = [y1 1-y1]; Psat = vp12(A,B,C,T); % data 
x10 = 0.7; % initial guess for mole fraction x 
[x1, fval] = fzero(@dewpf,x10,[],y,T,Psat); x = [x1 1-x1]; 
gamma = gamma12(x1,T); P = sum(x.*gamma.*Psat); 
fprintf('(2) Dew P: P = %g, x1 = %g, x2 = %g\n',P,x(1),x(2)); 
fprintf('gamma1 = %g, gamma2 = %g\n',gamma(1),gamma(2)); 
% (3) Bubble T calculation 
P = 101.33; x1 = 0.85; x = [x1 1-x1]; % data 
T0 = 300; % initial guess for temperature 
[T, fval] = fzero(@bubtf,T0,[],x1,A,B,C,P); 
gamma = gamma12(x1,T); Psat = vp12(A,B,C,T); y = (x.*gamma.*Psat)/P; 
fprintf('(3) Bubble T: T = %g, y1 = %g, y2 = %g\n',T,y(1),y(2)); 
fprintf('gamma1 = %g, gamma2 = %g\n',gamma(1),gamma(2)); 
% (4) Dew T calculation 
P = 101.33; y1 = 0.40; y = [y1 1-y1]; % data 
T0 = 330; % initial guess for temperature 
Psat0 = vp12(A,B,C,T0); x1 = y1*P/Psat0(1); crx = 1; cx = 1e-6; 
while crx > cx    
[T, fval] = fzero(@dewfun,T0,[],y1,x1,A,B,C,P);    
gamma = gamma12(x1,T); Psat = vp12(A,B,C,T);    
x = y*P./(gamma.*Psat); crx = abs(x - x1); x1 = x(1); 
end 
fprintf('(4) Dew T: T = %g, x1 = %g, x2 = %g\n',T,x(1),x(2)); 
fprintf('gamma1 = %g, gamma2 = %g\n',gamma(1),gamma(2));  

function gamma = gamma12(x1,T) % Calculation of activity coefficients  
A = 2.771 - 0.00523*T; % T(K) 
gamma(1) = exp(A*(1-x1)^2); gamma(2) = exp(A*x1^2); 
end 
function Psat = vp12(A,B,C,T) % Vapor pressure by Antoine equation (T: K) 
Psat = exp(A - B./(T - C));  
end 
function f = dewpf(x1,y,T,Psat) % Define equation for (2) dew P calculation 
x = [x1 1-x1]; gamma = gamma12(x1,T); P = sum(x.*gamma.*Psat); 
f = sum(y*P./(gamma.*Psat)) - 1; 
end 
function f = bubtf(T,x1,A,B,C,P) % Define equation for (3) Bubble T calculation 
x = [x1 1-x1]; gamma = gamma12(x1,T); Psat = vp12(A,B,C,T); 
f = sum(x.*gamma.*Psat) - P; 
end 
function fp = dewfun(T,y1,x1,A,B,C,P) % Define equation for (4) Dew T calculation 
y = [y1 1-y1]; Psat = vp12(A,B,C,T); gamma = gamma12(x1,T); 
fp = sum(y*P./(gamma.*Psat)) - 1; 
end 