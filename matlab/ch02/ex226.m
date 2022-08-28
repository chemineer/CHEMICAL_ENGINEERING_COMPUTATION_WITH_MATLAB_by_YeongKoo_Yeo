% nonreg.m: nonlinear regression 
clear all; 
t  =  [2.1 2.2 2.4 2.5 2.6 2.8 3.0 3.1]; z  =  [min(t): 0.01:max(t)]; 
v  =  [3.091 2.699 1.801 1.698 1.412 1.297 0.702 0.597]; 
X  =  [exp(t); t.*exp(t)]'; Y  =  [1./v]'; 
C  =  inv(X'*X)*X'*Y; A  =  1/C(1), B  =  A*C(2) % regression by least squares 
yz  =  A*exp(-z)./(1 + B*z); 
f  =  @(D,t) D(1)*exp(-t)./(1 + D(2)*t); % model to be used in the built-in function nlinfit 
D0  =  [1 1]'; % initial guess for parameter vector 
D  =  nlinfit(t,v,f,D0); % use of built-in function nlinfit 
plot(z,yz,z,f(D,z),'--',t,v,'o'), xlabel('t'), ylabel('v') 
legend('Nonlinear regression','Built-in fun nlinfit','Data','location', 'best')  